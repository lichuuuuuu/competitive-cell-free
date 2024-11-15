function [Hhat,H,B] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p)
%Generate the channel realizations and estimates of these channels for all
%UEs in the entire network. The channels are modeled as correlated
%Rayleigh fading and the MMSE estimator is used.
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
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k 
%                    in setup n, normalized by the noise power
%nbrOfRealizations = Number of channel realizations
%L                 = Number of APs
%K                 = Number of UEs in the network
%N                 = Number of antennas per AP
%tau_p             = Number of orthogonal pilots
%pilotIndex        = Vector containing the pilot assigned to each UE
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%Hhat         = Matrix with dimension L*N x nbrOfRealizations x K where
%               (:,n,k) is the estimated collective channel to UE k at
%               channel realization n.
%H            = Matrix with dimension L*N x nbrOfRealizations x K with the
%               true channel realizations. The matrix is organized in the
%               same way as Hhat_MMSE.
%B            = Matrix with dimension N x N x L x K where (:,:,l,j) is the
%               spatial correlation matrix of the estimate between AP l and
%               UE k in setup n, normalized by the noise power


%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
%各APのアンテナごとに無相関なレイリーフェージングチャネルHを生成
H = (randn(L*N,nbrOfRealizations,K)+1i*randn(L*N,nbrOfRealizations,K));


%Go through all channels and apply the spatial correlation matrices
%相関行列を適用することで任意の地理的または物理的条件下での信号伝播の特性をチャネルHに組み込む
%これはypの計算に使うだけ(受信したパイロット信号を再現するために使う)
% lはAPのループ    
for l = 1:L
    %APに関連するすべてのユーザー機器のループ        
    for k = 1:K
        
        %Apply correlation to the uncorrelated channel realizations
        %各APとUEのペアごとに定義された空間相関行列 R(:,:,l,k) の平方根を計算(この計算により、チャネルHに適用するための相関行列を作成)
        Rsqrt = sqrtm(R(:,:,l,k));
        %相関行列を掛けることで、無相関なチャネルHに所定の相関構造を導入
        H((l-1)*N+1:l*N,:,k) = sqrt(0.5)*Rsqrt*H((l-1)*N+1:l*N,:,k);
        
    end
    
end


%% Perform channel estimation

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));


%Prepare to store results
Hhat = zeros(L*N,nbrOfRealizations,K);

if nargout>2
    B = zeros(size(R));
end


%Go through all APs
%各APにおいて、パイロット信号に基づいてMMSE（最小平均二乗誤差）推定でチャネルHを推定する
% lはAPのループ  ．
for l = 1:L
    
    %Go through all pilots
    %tはパイロット信号のループ 
    for t = 1:tau_p
        
        %Compute processed pilot signal for all UEs that use pilot t
        %パイロット信号tを使用するすべてのUEからの受信信号の合計を計算
        %sqrt(p)*tau_p はUEからの送信電力とパイロット信号の長さを考慮したスケーリング係数
        %Np(:,:,l,t) はノイズの成分
        yp = sqrt(p)*tau_p*sum(H((l-1)*N+1:l*N,:,t==pilotIndex),3) + sqrt(tau_p)*Np(:,:,l,t);
        
        %Compute the matrix that is inverted in the MMSE estimator
        %MMSEの計算に必要な共分散行列の逆行列PsiInvを求める
        %ユーザー機器（UE）のAP送信電力，パイロット信号の長さ，パイロット信号tを使用するUEのAPlに対する空間相関行列の合計の数列と単位行列の和
        PsiInv = (p*tau_p*sum(R(:,:,l,t==pilotIndex),4) + eyeN);
        
        %Go through all UEs that use pilot t
        %kはUEのループ 
        %パイロットインデックスが t と等しい全てのUEのインデックスを返す
        for k = find(t==pilotIndex)'
            
            %Compute the MMSE estimate
            %AP l とUE k 間の空間相関行列 R(:,:,l,k) と　共分散行列の逆行列 PsiInvの割り算でウェイトRPsiを求める
            RPsi = R(:,:,l,k) / PsiInv;
            %処理されたパイロット信号ypをUEの送信電力sqrt(p)でスケーリングされたRPsiでフィルタリングし，チャネル行列H(AP lとUE k 間のチャネル)を推定
            Hhat((l-1)*N+1:l*N,:,k) = sqrt(p)*RPsi*yp;
            
            %Compute the spatial correlation matrix of the estimate
            %チャネルの空間相関行列Bを計算(MMSE推定の統計的特性) 
            if nargout>2
                B(:,:,l,k) = p*tau_p*RPsi*R(:,:,l,k);
            end
            
        end
        
    end
    
end
