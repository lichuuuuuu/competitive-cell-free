function [SE_MR,SE_MMSE,sumSE_SIC] = functionComputeSE_AP_uplink(Hhat,H,R,B,tau_c,tau_p,nbrOfRealizations,N,K,L,p,p1)
%Compute uplink SE for Cell-free mMIMO for the four different receiver
%cooperation levels, using either MR or MMSE/L-MMSE combining
%
%This function was developed as a part of the paper:
%
%Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
%Competitive With MMSE Processing and Centralized Implementation,"
%IEEE Transactions on Wireless Communications, To appear.
%
%Download article: https://arxiv.org/abs/1903.10611
%
%This is version 1.0 (Last edited: 2019-03-19)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%
%INPUT:
%Hhat              = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel from
%                    all BSs to UE k at channel realization n.
%H                 = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the true collective channel from all
%                    BSs to UE k at channel realization n.
%R                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix between BS
%                    l and UE k in setup n, normalized by the noise power
%B                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix of the
%                    estimate between BS l and UE k in setup n, normalized
%                    by the noise power
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences and number of UEs per cell
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Total number of UEs
%L                 = Number of APs
%p                 = Matrix K x 1 where element k is the uplink transmit
%                    power of UE k (If it is a scalar, the same value is
%                    used for all users)
%p1                = (Optional) Same as p but only for Level 1
%
%OUTPUT:
%SE_MR     = K x 4 matrix where the (k,n):th element is the uplink SE of 
%            UE k achieved with MR combining at cooperation level n
%SE_MMMSE  = Same as SE_MR but with MMSE or L-MMSE combining
%sumSE_SIC = Scalar with the sum SE achieved using MMSE-SIC combining



%If only one transmit power is provided, use the same for all the UEs
if length(p) == 1
   p = p*ones(K,1);
end

%If no specific Level 1 transmit powers are provided, use the same as for
%the other levels
if nargin<12
    p1 = p;
end


%Store identity matrices of different sizes
eyeN = eye(N);
eyeK = eye(K);
eyeNL = eye(N*L);

%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_MR = zeros(K,4);
SE_MMSE = zeros(K,4);

%Compute sum of all estimation error correlation matrices at every BS
%誤差共分散行列の合計を計算
%通常の送信パワー pを使用した際の、全UEに対する誤差共分散行列の合計を保存するための行列
C_tot = zeros(N,N,L);
%レベル1の送信パワー pを使用した際の、全UEに対する誤差共分散行列の合計を保存するための行列
C1_tot = zeros(N,N,L);
%kはUEのループ
for k = 1:K
    %UEkについて、実際の空間相関行列 R に各UEの送信パワー p(k) もしくは p1(k) を掛けた値と推定された空間相関行列 B の差を取り、UEkに由来する誤差の共分散の寄与を計算
    C_tot = C_tot + p(k)*(R(:,:,:,k)-B(:,:,:,k));
    C1_tot = C1_tot + p1(k)*(R(:,:,:,k)-B(:,:,:,k));
end

%システム全体の誤差共分散行列を格納するための大きなブロック対角行列
C_tot_blk = zeros(N*L,N*L);

for l = 1:L
    %システム全体の誤差共分散行列を構築する
    % 各APlに対して、C_tot(:,:,l) を C_tot_blk の対応するブロック位置に代入.
    C_tot_blk(1+(l-1)*N:l*N,1+(l-1)*N:l*N) = C_tot(:,:,l);
end

%Diagonal matrix with transmit powers and its square root
Dp = diag(p);
Dp12 = diag(sqrt(p));
Dp112 = diag(sqrt(p1));


%Prepare to save simulation results
SE_MR_level1 = zeros(K,L);
SE_MMSE_level1 = zeros(K,L);

signal_MR_level23 = zeros(L,K);
scaling_MR_level23 = zeros(L,K);
Gp_MR_level23 = zeros(L,L,K);

signal_MMSE_level23 = zeros(L,K);
scaling_MMSE_level23 = zeros(L,K);
G_MMSE_level23 = zeros(L,L,K);

sumSE_SIC = 0;


%% Go through all channel realizations
%nはチャネルHのループ
for n = 1:nbrOfRealizations
    
    
    %Level 4
    
    %Extract channel estimate realizations from all UEs to all APs
    %n番目のチャネル実現でのすべてのUEからの推定チャネルを抽出
    Hhatallj = reshape(Hhat(:,n,:),[N*L K]);
    
    %Compute MR combining
    %推定されたチャネル行列 Hhatallj をそのまま最大比結合の結合行列として使用
    V_MR = Hhatallj;
    
    %Compute MMSE combining
    %MMSE結合ベクトルを計算
    %チャネル推定行列の自己共分散(チャネル推定(Hhatallj)*送信電力の逆行列(Dp)*チャネル推定の逆行列(Hhatallj))，システム全体の誤差共分散行列，単位行列のノイズ(eyeNL)，の和で信号+ノイズ+干渉の全体的な共分散行列を導き，スケーリングされたチャネル推定行列で割る
    V_MMSE = ((Hhatallj*Dp*Hhatallj')+C_tot_blk+eyeNL)\(Hhatallj*Dp);
    

    %Go through all UEs
    %kはUEのループ
    for k = 1:K
        %MRを使って瞬時スペクトル効率を計算
        %%MR combining
        %UEkの最大比結合の結合行列を抽出(特定のユーザー機器（UE）から送信される信号を特定し、それを他のノイズや干渉信号から分離するために使用される)
        v = V_MR(:,k); %Extract combining vector
        
        %Compute numerator and denominator of instantaneous SINR at Level 4
        %瞬時SINR（信号対干渉ノイズ比）の計算
        %分子として．UE kに対する結合ベクトル v とそのUEのチャネル推定 Hhatallj(:,k) の内積の絶対値の二乗を計算
        numerator = p(k)*abs(v'*Hhatallj(:,k))^2;
        %分母として，全UEからの信号の合計パワーを計算し（norm(v'*Hhatallj*Dp12)^2）、そこから自UEの信号パワーを差し引き、ノイズと干渉の合計を計算する(v'*(C_tot_blk+eyeNL)*v は、推定誤差とノイズによる干渉の寄与)
        denominator = norm(v'*Hhatallj*Dp12)^2 + v'*(C_tot_blk+eyeNL)*v - numerator;
        
        %Compute instantaneous SE for one channel realization
        %瞬時スペクトル効率を計算する．上記で計算した瞬時SINRを使ってスペクトル効率を計算し，全試行数nbrOfRealizationsによって平均化
        %prelogFactorは、シンボルレートとピロットシンボルの使用によって調整された前置係数
        SE_MR(k,4) = SE_MR(k,4) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
        %MMSEを使って瞬時スペクトル効率を計算
        %%MMSE combining
        v = V_MMSE(:,k); %Extract combining vector
        
        %Compute numerator and denominator of instantaneous SINR at Level 4
        numerator = p(k)*abs(v'*Hhatallj(:,k))^2;
        denominator = norm(v'*Hhatallj*Dp12)^2 + v'*(C_tot_blk+eyeNL)*v - numerator;
        
        %Compute instantaneous SE for one channel realization
        SE_MMSE(k,4) = SE_MMSE(k,4) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        
    end
    
    
    %Compute sum SE with MMSE-SIC combining at Level 4
    %最小平均二乗誤差-逐次干渉キャンセリング（MMSE-SIC）結合を用いた場合の全体のスペクトル効率（SE）を計算
    if nargout > 2
        %(det(eyeK +
        %Dp*V_MR'*((C_tot_blk+eyeNL)\V_MR)))はシステム全体のSINR行列の行列式を計算，V_MR' はMR結合行列の共役転置，C_tot_blk + eyeNL は全ユーザーにわたるノイズと誤差共分散行列に単位行列を加えたもの
        sumSE_SIC = sumSE_SIC + prelogFactor*real(log2(det(eyeK + Dp*V_MR'*((C_tot_blk+eyeNL)\V_MR))))/nbrOfRealizations;
        
    end
    
    
    
    
    %Levels 1-3
    gp_MR_level23 = zeros(L,K,K);
    gp_MMSE_level23 = zeros(L,K,K);
    
    
    %Go through all APs
    for l = 1:L
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H(1+(l-1)*N:l*N,n,:),[N K]);
        
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(1+(l-1)*N:l*N,n,:),[N K]);
        
        
        %Compute MR combining
        V_MR = Hhatallj;
        V_MMSE = ((Hhatallj*Dp*Hhatallj')+C_tot(:,:,l)+eyeN)\(V_MR*Dp);
        

        %Go through all UEs
        for k = 1:K
            
            %MRを使って瞬時スペクトル効率を計算
            %%MR combining
            v = V_MR(:,k); %Extract combining vector
            
            
            %Level 2 and Level 3
            %結合ベクトル v とUE k からAP l までの実際のチャネル Hallj の内積を取り、特定のUEからの信号成分を抽出し,平均化処理
            signal_MR_level23(l,k) = signal_MR_level23(l,k) + (v'*Hallj(:,k))/nbrOfRealizations;
            %結合ベクトル v と全UEからAP l までのチャネル Hallj の内積を取り、全UEの信号を同時に抽出
            gp_MR_level23(l,:,k) = gp_MR_level23(l,:,k) + (v'*Hallj)*Dp12;
            %結合ベクトル v のノルムの二乗を計算し，チャネル実現の数で平均化．信号処理におけるエネルギーの大きさを求める
            scaling_MR_level23(l,k) = scaling_MR_level23(l,k) + norm(v).^2/nbrOfRealizations;
            
            
            %Level 1
            
            %Compute numerator and denominator of instantaneous SINR
            numerator = p1(k)*abs(v'*Hhatallj(:,k))^2;
            denominator = norm((v'*Hhatallj)*Dp112)^2 + v'*(C1_tot(:,:,l)+eyeN)*v - numerator;
            
            %Compute instantaneous SE for one channel realization
            SE_MR_level1(k,l) = SE_MR_level1(k,l) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
            
            
            %MMSEを使って瞬時スペクトル効率を計算
            %%MMSE combining
            v = V_MMSE(:,k); %Extract combining vector
            
            
            %Level 2 and Level 3
            %UEに対するAPからの信号の合計応答が計算される(結合ベクトル v と実際のチャネル Hallj との内積)
            signal_MMSE_level23(l,k) = signal_MMSE_level23(l,k) + (v'*Hallj(:,k))/nbrOfRealizations;
            %グローバルパイロット応答の集約が計算される(結合ベクトル v と実際のチャネル Hallj との内積に各UEの送信パワーの平方根を表す Dp12（対角行列）を乗算し，送信パワーに応じた信号応答の調整を行う)
            gp_MMSE_level23(l,:,k) = gp_MMSE_level23(l,:,k) + (v'*Hallj)*Dp12;
            %スケーリングパラメータの集約が計算される．結合ベクトル v のノルムの二乗を計算(結合ベクトルがどの程度のエネルギーを信号に適用しているかを示す)
            scaling_MMSE_level23(l,k) = scaling_MMSE_level23(l,k) + norm(v).^2/nbrOfRealizations;
            
            
            %Level 1
            
            %Compute numerator and denominator of instantaneous SINR
            numerator = p1(k)*abs(v'*Hhatallj(:,k))^2;
            denominator = norm((v'*Hhatallj)*Dp112)^2 + v'*(C1_tot(:,:,l)+eyeN)*v - numerator;
            
            %Compute instantaneous SE for one channel realization
            %特定のチャネル実現におけるユーザー機器（UE）のスペクトル効率（SE）を計算し、それを集約
            SE_MMSE_level1(k,l) = SE_MMSE_level1(k,l) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
            
            
        end
        
    end
    
    %Compute averaging of terms for Level 2 and Level 3
    for k = 1:K
        %通信チャネル上のすべての実現にわたって収集されたデータを用いて、干渉の統計的特性を評価
        %gp_MR_level23(:,:,k)は最大比結合（MR）で得られる信号の結合成果を格納しているもの，これの自己相関行列を求めて信号のエネルギー分布と干渉の統計的特性を明らかにする
        Gp_MR_level23(:,:,k) = Gp_MR_level23(:,:,k) + gp_MR_level23(:,:,k)*gp_MR_level23(:,:,k)'/nbrOfRealizations;
        %gp_MMSE_level23(:,:,k)は最小平均二乗誤差結合（MMSE）で得られる信号の結合成果を格納しているもの，これの自己相関行列を求めて信号のエネルギー分布と干渉の統計的特性を明らかにする
        G_MMSE_level23(:,:,k) = G_MMSE_level23(:,:,k) + gp_MMSE_level23(:,:,k)*gp_MMSE_level23(:,:,k)'/nbrOfRealizations;
        
    end
    
end


%Compute SE for Level 1
%集約したMRとMMSEのペクトル効率（SE）のうち最大のものをスペクトル効率（SE）として抽出
SE_MR(:,1) = max(SE_MR_level1,[],2);
SE_MMSE(:,1) = max(SE_MMSE_level1,[],2);


%Compute SE for Level 2 and Level 3
%MRとMMSEの瞬時スペクトル効率を最終的に計算
for k = 1:K
    
    %With MR combining
    
    b = signal_MR_level23(:,k);
    %各UEの結合信号の共分散行列
    A = Gp_MR_level23(:,:,k) + diag(scaling_MR_level23(:,k)) - p(k)*(b*b');
    %逆行列計算（A\b）を用いて、Level 3でのMRのスペクトル効率を計算
    SE_MR(k,3) = prelogFactor*real(log2(1+p(k)*b'*(A\b))); 
    %全APからの信号の平均強度を用いて、より単純な計算（mean(b)）でLevel 2でのMRのスペクトル効率を計算
    SE_MR(k,2) = prelogFactor*real(log2(1+p(k)*abs(mean(b)).^2 / mean(mean(A))));   

    %With L-MMSE combining

    b = signal_MMSE_level23(:,k);
    %各UEの結合信号の共分散行列
    A = G_MMSE_level23(:,:,k) + diag(scaling_MMSE_level23(:,k)) - p(k)*(b*b');
    %逆行列計算（A\b）を用いて、Level 3でのMMSEのスペクトル効率を計算
    SE_MMSE(k,3) = prelogFactor*real(log2(1+p(k)*b'*(A\b)));  
    %全APからの信号の平均強度を用いて、より単純な計算（mean(b)）でLevel 2でのMMSEのスペクトル効率を計算
    SE_MMSE(k,2) = prelogFactor*real(log2(1+p(k)*abs(mean(b)).^2 / mean(mean(A))));   

end
