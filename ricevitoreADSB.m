%Tesi - Laurea Triennale ingegneria Fisica
%
%Titolo: Decoding Live Wireless Signals with MATLAB and RTL-SDR
%
%Scritto da Lorenzo Maria Scarrone
%
%Corso di riferimento a scelta
%         INTERNET E COMUNICAZIONI
%
%referente: prof. Andrea Carena


%CARICO LA CONFIGURAZIONE PER IL RICEVITORE ADSB
adsbParam = ConfigAdsbReciverAndStandardParam();
dataCursor = 1;

%SCANSIONE PACCHETTI ADSB (ETERE)                  
%ADSBrx = RTLreciver(adsbParam);

%CARICO PACCHETTO DATI PREACQUISITO (offline)
load('acquisizionisalvate/rxADSB_3aerei.mat');

for cur=1:adsbParam.NBlocks
    %sovracampionamento del segnale ricevuto: ottengo SpS intero
    % 2.4MHz x 5 = 12MHz --> 12MHz/(1Mb/s) = 12 SpS
    rxFil = adsbInterpolator(ADSBrx(cur,:)',adsbParam);
    %normalizzo il segnale interpolato   
    rxAbs = abs(rxFil)/max(abs(rxFil));
    
    %ottengo i pacchetti dati sincronizzando il flusso dati con il preamble    
    [adsbPacket,errors] = CorrelateAndSync( rxAbs ,adsbParam);
    
    %ANALISI DEI PACCHETTI ESTRATTI
    for ii=1:size(adsbPacket,1)
        %VERIFICO L'ESISTENZA DEL PACCHETO
        if errors(ii,1) ~= 1
            %SI PROCEDE CON LA DEMODULAZIONE DEL PACCHETTO
            adsbBits= adsbDemod(adsbPacket(ii,:)',adsbParam);
            %VERIFICA DELLA VALIDITA' DEL PACCHETTO
            invalidPacket = CRCcheck(adsbBits);
            if invalidPacket == false
                %SUDDIVISIONE DEL PACCHETTO PER LA DECODIFICA
                ADSBdata(dataCursor) = decodeData(adsbBits);
                dataCursor = dataCursor+1;
            end
        end
    end        
end
%DECODIFICA DEI PACCHETTI ESTRATTI
if exist('ADSBdata')
    DecodedData = decodeHexData(ADSBdata);
    struct2table(DecodedData)
    plotData(DecodedData);  
else
    warning('Nessun pacchetto ADS-B rilevato')
end

%-----------------------------. . . . . . . . . ---------------------------
%-------------------------- . . .  FUNZIONI . . . -------------------------
%-----------------------------. . . . . . . . . ---------------------------

function plotData(decData)
    cellData = struct2cell(decData);
    ICAO_list = unique(cellData(1,1,:));

    figure();
    for i=1:length(ICAO_list)
        n=1;
        for j=1:length(cellData)
            if isequal(ICAO_list(i),cellData(1,1,j))
                x(n) = cellData(2,1,j);
                y(n) = cellData(3,1,j);
                h(n) = cellData(4,1,j);
                n=n+1;
            end
        end
        x1=cell2mat(x);y1=cell2mat(y);
        plot(x1,y1,'-o');
        hold on
        plot(x1,y1);
        for k=1:length(x1)
            txt = "icao: "+cell2mat(ICAO_list(1))+" altitude: "+num2str(cell2mat(h(k)));
            text(x1(k),y1(k),txt,'HorizontalAlignment','right');
        end
        hold on
        clear x y
    end
end

%FUNZIONE PER ESTRARRE I PACCHETTI DA FRAME ACQUISITO
function [sync,error] = CorrelateAndSync(Signal, param)    
    SyncSeq = param.preamble; %CARICO LA SYNC SEQUENCE
    SpS = param.SpS; 
    %CALCOLO LA CORRELAZIONE PER OTTENERE I PICCHI DA CUI ESTRARRE I DATI
    [Xcor,Lag] = xcorr(SyncSeq,Signal);
    %RICERCA DEI PICCHI DI CORRELAZIONE CON IL TEST DI SOGLIA
    [~,dl] = findpeaks(Xcor/max(Xcor),Lag,'MinPeakHeight',0.8);
    %SALVO LA POSIZIONE DEI PICCHI NEL FRAME
    Pos = abs(dl);
    %RICERCA DEI PACCHETTI
    for ii=1:length(Pos)
        start_pos = Pos(ii)+length(SyncSeq)+1;
        error(ii,1)=0;
        %VERIFICO CHE LA DIMENSIONE DEL PACCHETTO NON ESCA DAL FRAME
        if length(Signal) >= start_pos+ param.LongPacketNumBits*SpS
            %RESTITUISCO IL PACCHETTO NUMERATO
            sync(ii,:) = abs(Signal(start_pos:start_pos+param.LongPacketNumBits*SpS-1));
        else
            error(ii,1)=1;
            sync(ii,:) = (1:1344)*0;
        end
    end
    
end

%FUNZIONE INTERPOLATRICE PER AUMENTARE IL NUMERO DI CAMPIONI
function z = adsbInterpolator(y, param)
%rendo la funzione InterpFil persistente = non deve essere rigenerata tutte
%le volte
    persistent interpFil
    if isempty(interpFil)
        %USO FILTRO FIR INTERPOLATOR PER AUMENTARE IL NUMERO DI CAMPIONI
        %ref; https://it.wikipedia.org/wiki/Finite_impulse_response
        interpFil = dsp.FIRInterpolator(param.InterpolationFactor,...
            param.InterpolationFilterCoefficients);
    end
    if param.InterpolationFactor > 1
        z = interpFil(y);
    else
        z=y;
    end
end

%FUNZIONE PER LA DECODIFICA DEL PACCHETTO
function [structData,error] = decodeData(BitVector)
    bits = logical(BitVector)';
    %Pacchetto ADSB vettore esadecimale completo
    structData.hexVec = binaryVectorToHex(bits);
    
    %Composizione del pacchetto ADSB Extended Squitter
    %[DF 1:5][Capability 6:8][ICAO 9:32][DATA 33:88][CRC 89:112]
    %Composizione del pacchetto ADSB Short Squitter
    %[DF 1:5][Capability 6:8][ICAO 9:32][CRC 33:56]
    
    structData.DF = bi2de(bits(1:5),'left-msb');
    structData.CA = bi2de(bits(6:8),'left-msb');
    structData.ICAO = binaryVectorToHex(bits(9:32));
    
    if structData.DF == 17 %Extended squitter
        structData.DATA = binaryVectorToHex(bits(33:88));
        structData.TC = bi2de(bits(33:37),'left-msb');
        structData.CPRf = bits(54); %cpr flag even odd
        structData.PI = binaryVectorToHex(bits(89:112));
    elseif structData.DF == 11 %Short squitter
        structData.PI = binaryVectorToHex(bits(33:56));
    else
        error = 1;
    end   
end

%FUNZIONI PER LA DECODIFICA DEI DATI DEL PACCHETTO
function structData = decodeHexData(ADSBdataStruct)
%CARICO TUTTI I FLAG CPRf e TC DEI PACCHETTI OTTENUTI DA decodeData
    CPRf = [ADSBdataStruct.CPRf];
    TC   = [ADSBdataStruct.TC];
    cursor = 1;
    for ii=1:length(TC)-1
        for jj=ii+1:length(TC)
%PER LA DECODIFICA DELLA POSIZIONE E' ESSENZIALE CHE I PACCHETTI DATI
%APPARTENGANO ALLO STESSO AEREO QUINDI RICHIAMO IN MEMORIA ICAO E ICAO1
%RIFERITI AL PRIMO E SECONDO PACCHETTO ANALIZZATO E IL VETTORE ADSB PER
%ESTRARRE I DATI
            ICAO = ADSBdataStruct(ii).ICAO;
            hBITS = ADSBdataStruct(ii).hexVec;           
            ICAO1 = ADSBdataStruct(jj).ICAO;
            hBITS1 = ADSBdataStruct(jj).hexVec;
            %VERIFICA ID AEREO
            if (ICAO == ICAO1)
                %VERIFICA DEL CPR FLAG E DEL TC PER OTTENERE LA POSIZIONE
                if (CPRf(ii)==0)&&(CPRf(jj) == 1)&&(TC(ii)==TC(jj))&&(TC(ii)>=9 && TC(ii)<=18)
                    
                    BITS = hexToBinaryVector(hBITS);
                    BITS1 = hexToBinaryVector(hBITS1);
                    
                    [Latitude,Longitude,Altitude] = getADSBvalue(BITS,BITS1);
                    %salvo i dati ottenuti
                    structData(cursor).ICAO = ICAO;
                    structData(cursor).Latitude = Latitude;
                    structData(cursor).Longitude = Longitude;
                    structData(cursor).Altitude = Altitude;
                    
                    cursor = cursor+1;                    
                end
            end
        end
    end
end

function [latitude,longitude,altitude] = getADSBvalue(BITS,BITS1)
                NZ = 15;
                dLat_even = 360/(4*NZ);
                dLat_odd = 360/(4*NZ-1);
                %acquisizione dei dati da valutare
                cpr_lat_even = bi2de(BITS(55:71),'left-msb')/131072;
                cpr_lon_even = bi2de(BITS(72:88),'left-msb')/131072;
                cpr_lat_odd = bi2de(BITS1(55:71),'left-msb')/131072;
                cpr_lon_odd = bi2de(BITS1(72:88),'left-msb')/131072;  
                
                %calcolo la latitudine
                j = floor(59*cpr_lat_even-60*cpr_lat_odd+0.5);
                
                Lat_even = dLat_even*(mod(j,60)+cpr_lat_even);
                Lat_odd = dLat_odd*(mod(j,59)+cpr_lat_odd);
                
                if Lat_even >= 270
                    Lat_even = Lat_even - 360;
                end
                if Lat_odd >= 270
                    Lat_odd = Lat_odd - 360;
                end                
                latitude = Lat_odd;
                
                %calcolo la longitudine
                ni = max(NL(cpr_lat_odd)-1,1);
                dLon = 360/ni;
                
                m = floor(cpr_lon_even*(NL(Lat_odd)-1)-cpr_lon_odd*NL(Lat_odd)+0.5);
                Lon_odd = dLon*(mod(m,ni)+cpr_lon_odd);
                
                if (Lon_odd>=180)
                    Lon_odd = Lon-360;
                end                
                longitude = Lon_odd;
                
                %calcolo l'altitudine
                altitude = bi2de([BITS1(41:47) BITS1(49:52)],'left-msb');
                if BITS1(48)==1
                    q_bit_mult=25;
                else
                    q_bit_mult=100;
                end;
                altitude = altitude*q_bit_mult-1000;
end

function nl = NL(latitude)
    NZ = 15;
    nl = floor(2*pi/acos(1-(1-cos(pi/(2*NZ)/cos(pi/180*latitude)^2))));
end

%FUNZIONE RICEVITORE E ACQUISIZIONE DATI                 
function rxSig = RTLreciver(RTLParam)

    rxADSb = comm.SDRRTLReceiver('0',...
        'CenterFrequency',RTLParam.CenterFrequency,...
        'EnableTunerAGC',false,...
        'TunerGain',RTLParam.TunerGain,...
        'SampleRate',RTLParam.SampleRate,...
        'OutputDataType','single',...
        'SamplesPerFrame',RTLParam.NsBlock,... 
        'FrequencyCorrection',0);

    if ~isempty(sdrinfo(rxADSb.RadioAddress))
          rxSig = zeros(RTLParam.NBlocks,RTLParam.NsBlock);
        for counter=1:RTLParam.NBlocks        
            rxSig(counter,:)=step(rxADSb)';
        end
    else
        warning('No disp sdr rilevato')
    end

    release(rxADSb);
end

%DEMODULAZIONE DEL SEGNALE - PASSO IN BINARIO
function z = adsbDemod(y,param)
%LA SPIEGAZIONE DEL FUNZIONAMENTO DI QUESTA FUNZIONE E CONTENUTA NEL REPORT
%IN MODO CHIARO
%differenza tra le due modulanti (1,0)-(0,1) = (1,-1)
    diffS_12 = [ones(param.SpC,1); -ones(param.SpC,1)];
%ottengo il la matrice [numBits x SpS] ovvero [112 x 12]
    numBits = size(y,1) / param.SpS;
    yTemp = reshape(y, param.SpS, numBits)';
%prodotto tra [112 x 12 ] [12 x 1] = [112 x 1] vettore di bit dove 1 logico
%è maggiore di 0 mentre lo 0 logico è minore di 0
    ySoft = yTemp*diffS_12;
    z = uint8(ySoft > 0);    
end

%CONTROLLO CRC PACCHETTI PER LA LORO VALIDITA'
function error = CRCcheck(data)
    
    generator = logical([1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0 0 1 0 0 1]);
    msg = logical(data)';
    
    for cursor=1:88 %112 bit -24bit per la parità
        if logical(msg(cursor)) == 1
            msg(cursor:cursor+24) = xor(msg(cursor:cursor+24),generator);
        end
    end
    CRC = msg(89:end);
    
    if CRC == false
        error = false;
    else
        error = true;
    end
end

%ADS-B STANDARD PARAM
function adsbParam = ConfigAdsbReciverAndStandardParam()
% definizione dei parametri standard per un ricevitore adsb fatto con
% rtl-sdr
    %PARAMETRI RADIO
    adsbParam.CenterFrequency = 1090e6;
    adsbParam.SampleRate = 2.4e6;
    adsbParam.NsBlock = 51840;
    adsbParam.TunerGain = 60;
    adsbParam.NBlocks = 1000;
    %PARAMETRI SIMBOLI
    adsbParam.NbitPerSymbol = 2;
    adsbParam.SpS = 12;
    adsbParam.SpC = 6; 
    %PARAMETRI PACCHETTI
    adsbParam.LongPacketNumBits = 112;
    adsbParam.ShortPacketNumBits = 56;
    %SEQUENZA DI SYNC
    adsbParam.preamble = [ones(1,6) zeros(1,6) ones(1,6) zeros(1,6) ...
                          zeros(1,18) ones(1,6) zeros(1,6) ...
                          ones(1,6) zeros(1,36)]; 

    % PRESA DA ADSBExample/helperAdsbConfig.m                  
    adsbParam.InterpolationFilterCoefficients = ...
                    single(rcosdesign(0.5, 3, double(adsbParam.SpC)));
    % IL FATTORE DI INTERPOLAZIONE E' 5 PER PASSARE DA 2.4 A 12 SpS
    adsbParam.InterpolationFactor=5;
end