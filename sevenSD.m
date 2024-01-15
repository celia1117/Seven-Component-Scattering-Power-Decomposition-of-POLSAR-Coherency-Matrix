function [Ps,Pd,Pv,Ph,Pod,Pcd,Pmd] = sevenSD(T)
% 7成分目标分解

[lines,pixels,~,~] = size(T);
t11 = T(:,:,1,1); t12 = T(:,:,1,2); t13 = T(:,:,1,3);
t21 = T(:,:,2,1); t22 = T(:,:,2,2); t23 = T(:,:,2,3);
t31 = T(:,:,3,1); t32 = T(:,:,3,2); t33 = T(:,:,3,3);

TP = t11+t22+t33;


%% 计算初始的各散射分量
Ph = 2*abs(imag(t23));      % 螺旋散射
Pod = 2*abs(real(t13));     % 定向偶极子散射
Pmd = 2*abs(real(t23));     % 混合偶极子散射
Pcd = 2*abs(imag(t13));     % 耦合偶极子散射
Pv = 2*t33-Ph-Pmd-Pod-Pcd;  % 基本体散射

%% 计算条件C1
C1 = t11-t22-(7/8)*t33+(1/16)*(Pmd+Ph)-(15/16)*(Pcd+Pod);
Ph(C1>0) = NaN;
Pv(C1>0) = NaN;

%% 计算不同分支条件下的各散射模型分量
Ps = zeros(lines,pixels); Pd = zeros(lines,pixels);
Pmw = zeros(lines,pixels); Pcw = zeros(lines,pixels);
S = zeros(lines,pixels); D = zeros(lines,pixels);
C = zeros(lines,pixels);

for i = 1 : lines
    for j = 1 : pixels
        if C1(i,j)>0
            % 判断Pv是否大于0
            if (Pv(i,j)>0)
            else
                Ph(i,j) = 0;
                Pod(i,j) = 0;
                Pcd(i,j) = 0;
                Pmd(i,j) = 0;
                Pv(i,j) = 2*t33(i,j);
            end
            % 计算条件ratio
            ratio = 10*log10((t11(i,j)+t22(i,j)-2*real(t12(i,j)))/(t11(i,j)+t22(i,j)+2*real(t12(i,j))));
            
            % 根据ratio选择不同的体散射模型
            if (ratio<=-2)
                Pv(i,j) = (15/8)*Pv(i,j);
                S(i,j) = t11(i,j)-Pv(i,j)/2-Pod(i,j)/2-Pcd(i,j)/2;
                D(i,j) = t22(i,j)-(7/30)*Pv(i,j)-Ph(i,j)/2-Pmd(i,j)/2;
                C(i,j) = t12(i,j)-Pv(i,j)/6;
            elseif (ratio>-2)&&(ratio<2)
                Pv(i,j) = (2)*Pv(i,j);
                S(i,j) = t11(i,j)-Pv(i,j)/2-Pod(i,j)/2-Pcd(i,j)/2;
                D(i,j) = t22(i,j)-(1/4)*Pv(i,j)-Ph(i,j)/2-Pmd(i,j)/2;
                C(i,j) = t12(i,j);
            elseif (ratio>2)
                Pv(i,j) = (15/8)*Pv(i,j);
                S(i,j) = t11(i,j)-Pv(i,j)/2-Pod(i,j)/2-Pcd(i,j)/2;
                D(i,j) = t22(i,j)-(7/30)*Pv(i,j)-Ph(i,j)/2-Pmd(i,j)/2;
                C(i,j) = t12(i,j)+Pv(i,j)/6;
            end
            
            % 判断是否超过总能量
            if (Pv(i,j)+Pmd(i,j)+Ph(i,j)+Pcd(i,j)+Pod(i,j))>TP(i,j)
                Ps(i,j) = 0;
                Pd(i,j) = 0;
                Pv(i,j) = TP(i,j)-Ph(i,j)-Pmd(i,j)-Pod(i,j)-Pcd(i,j);
            else
                % 计算条件C0
                C0 = 2*t11(i,j)+Ph(i,j)+Pmd(i,j)-TP(i,j);
                % 根据C0判断
                if (C0>0)
                    Ps(i,j) = S(i,j) + (abs(C(i,j))^2)/S(i,j);
                    Pd(i,j) = D(i,j) - (abs(C(i,j))^2)/S(i,j);
                else
                    Pd(i,j) = S(i,j) - (abs(C(i,j))^2)/D(i,j);
                    Ps(i,j) = D(i,j) + (abs(C(i,j))^2)/D(i,j);
                end
                
                % 判断能量是否小于0
                if (Ps(i,j)>0)&&(Pd(i,j)<0)
                    Pd(i,j) = 0;
                    Ps(i,j) = TP(i,j)-Pv(i,j)-Ph(i,j)-Pmd(i,j)-Pod(i,j)-Pcd(i,j);
                elseif (Pd(i,j)>0)&&(Ps(i,j)<0)
                    Ps(i,j) = 0;
                    Pd(i,j) = TP(i,j)-Pv(i,j)-Ph(i,j)-Pmd(i,j)-Pod(i,j)-Pcd(i,j);
                else
                end
            end
        else
            % C1<=0
            % 判断Pv是否大于0
            if (Pv(i,j)>0)
            else
                Pv(i,j) = 0;
                Pmw(i,j) = Pmd(i,j) + Ph(i,j);
                Pcw(i,j) = Pcd(i,j) + Pod(i,j);
                % 判断Pmw与Pcw两者的关系
                if (Pmw(i,j)>Pcw(i,j))
                    Pcw(i,j) = 2*t33(i,j)-Pmw(i,j);
                    if (Pcw(i,j)<0)
                        Pcw(i,j) = 0;
                    else
                    end
                    % 判断Pod与Pcd之间的关系
                    if (Pod(i,j)>Pcd(i,j))
                        Pcd(i,j) = Pcw(i,j)-Pod(i,j);
                        if (Pcd(i,j)<0)
                            Pcd(i,j) = 0;
                        else
                        end
                    else
                        Pod(i,j) = Pcw(i,j)-Pcd(i,j);
                        if (Pod(i,j)<0)
                            Pod(i,j) = 0;
                        else
                        end
                    end
                else
                    Pmw(i,j) = 2*t33(i,j)-Pcw(i,j);
                    if (Pmw(i,j)<0)
                        Pmw(i,j) = 0;
                    else
                    end
                    % 判断Pmd与Ph之间的关系
                    if (Pmd(i,j)>Ph(i,j))
                        Ph(i,j) = Pmw(i,j)-Pmd(i,j);
                        if (Ph(i,j)<0)
                            Ph(i,j) = 0;
                        else
                        end
                    else
                        Pmd(i,j) = Pmw(i,j)-Ph(i,j);
                        if (Pmd(i,j)<0)
                            Pmd(i,j) = 0;
                        else
                        end
                    end
                end
            end
            
            % 计算各散射分量
            Pv(i,j) = (15/16)*Pv(i,j);
            S(i,j) = t11(i,j)-Pod(i,j)/2-Pcd(i,j)/2;
            D(i,j) = t22(i,j)-(7/15)*Pv(i,j)-Ph(i,j)/2-Pmd(i,j)/2;
            C(i,j) = t12(i,j);
            
            Ps(i,j) = S(i,j)-(abs(C(i,j)))^2/D(i,j);
            Pd(i,j) = D(i,j)+(abs(C(i,j)))^2/D(i,j);
            
            % 判断能量是否小于0
            if (Ps(i,j)>0)&&(Pd(i,j)<0)
                Pd(i,j) = 0;
                Ps(i,j) = TP(i,j)-Pv(i,j)-Ph(i,j)-Pmd(i,j)-Pod(i,j)-Pcd(i,j);
            elseif (Pd(i,j)>0)&&(Ps(i,j)<0)
                Ps(i,j) = 0;
                Pd(i,j) = TP(i,j)-Pv(i,j)-Ph(i,j)-Pmd(i,j)-Pod(i,j)-Pcd(i,j);
            else
            end
                        
        end
    end
end

% 将所有负能量置零
Ps(Ps<0) = 0;
Pd(Pd<0) = 0;
Pv(Pv<0) = 0;
Ph(Ph<0) = 0;
Pod(Pod<0) = 0;
Pcd(Pcd<0) = 0;
Pmd(Pmd<0) = 0;


end

