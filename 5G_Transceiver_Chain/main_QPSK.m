clc;
tb_length= input("enter the transport block length=\n");
message=randi([0,1],1,tb_length); %message is generated
%-----------------CRC Generation--------------------%
CRC= [24 23 6 5 1 0];
generator = zeros(1, length(CRC)+1);
for i=1:length(CRC)
    generator(1,CRC(1,i)+1)=1;
end
generator=flip(generator); 
%----------generator polynomial is generated------------%

txed_data=append_crc(message,generator);%CRC is appended to Transport block
C=length(txed_data)/(8448-24);            %Kcb for BG1=8448
C=ceil(C); %Number of Segmented CodeBlocks
L=24;   %Number of CRC bits to be appended as per Standards
B1=length(txed_data) + C*L;   %B1=Effective payload length
K1=B1/C;  %expected code block size, which is actually incorrect 


%-----------------Zc value decision takes place-------------%
ZcVec = [2:16 18:2:32 36:4:64 72:8:128 144:16:256 288:32:384];
for i=1:length(ZcVec)
    if 22*ZcVec(i)>=K1
        Zc=ZcVec(i); %Zc value is found
        break;
    end
end
%-----------------appropriate Zc is found--------------%
K=22*Zc; %Correct code Block Size
no_of_filler=K-K1; % Number of filler Bits is calculated    
p=K1-L;
W=reshape(txed_data, [p,C]); %Segmentation is done & stored in Matrix W
W=transpose(W);
V=zeros(C,p+L);
%-----------------Code Block CRC is added-----------------%
for i=1:C
    V(i,:)=append_crc(W(i,: ),generator);
end
%CodeBlock CRC is added to each Segmented block and stored in 2-D Matrix V

F=-1*ones(C,no_of_filler);
%-----------------Filler Bits are added-------------------%
V=[V F];   %V is 2-D matrix with CB-CRC along with Filler Bits
V=transpose(V);

G=100*162*2;   %Total number of bits which can be transmitted for 100PRBs 
E=G/C;    %Number of bits allowed to be sent for each segmented code block
nbg = 1; % Base graph 1
nldpcdecits = 25; % Decode with maximum no of iteration
  
%-----------------LDPC encoding---------------------------%
noise_power = (10^-5);
ldpc_encoded = double(LDPCEncode(V,nbg));
%-----------------LDPC encoding done-----------------------%

rate_mat=zeros(E,C);

%-----------------Rate Matching--------------------------------%
for m=1:E
    rate_mat(m,:)=ldpc_encoded(m,:);
end
%-----------------Rate Matching done--------------------------------%

%-----------------Interleaving--------------------------------%
rate_match=transpose(rate_mat);
order=2;
a1=E/order;
A=zeros(a1,order);
B=zeros(C,E);
for i=1:C
    AB=rate_match(i,:);
    A=reshape(AB, [a1,order]);
    A=transpose(A);
    for h=1:order
        B(i,h:order:end)=A(h,:);
    end
end
%-----------------Interleaving done--------------------------------%

Interleaved = reshape(B', 1, []); % Concatenation

%---------------------Scrambling--------------------------------%
x1=zeros(1,31);
x1(1)=1;
x2=zeros(1,31);
x2(1:8)=1;
x3=zeros(1,G+1600);
for k=1:(1600+G)
    x3(k)=double(xor(x1(1),x2(1)));
    r=double(xor(x1(1),x1(4)));
    x1(1:30)=x1(2:31);
    x1(31)=r;
    s=double(xor(double(xor(double(xor(x2(1),x2(2))),x2(3))),x2(4)));
    x2(1:30)=x2(2:31);
    x2(31)=s;
end
x4(1:G)=x3(1601:end);
scrambled=double(xor(Interleaved,x4));
%---------------------Scrambling Done---------------------------%

%-----------------QPSK Modulation-------------%
Rate_Matched_new=-2*(scrambled-0.5);
mod_output=(Rate_Matched_new(1:2:end)+1j*Rate_Matched_new(2:2:end))/sqrt(2);
% mod_output=(Rate_Matched_new(1:2:end)+1j*Rate_Matched_new(2:2:end));
noise = sqrt(noise_power)*randn(size(mod_output));
inphase_quad_noise=(noise+1j*noise);
rx_sig = mod_output + inphase_quad_noise;
%--------------QPSK Modulation Done-----------%

%--------------QPSK De-Modulation-------------%

real_rx_sig=real(rx_sig);
img_rx_sig=imag(rx_sig);

llr0 =  abs(1/sqrt(2) + real_rx_sig);  
llr1 =  abs(-1/sqrt(2) + real_rx_sig);    
llr_inphase = log(llr0./llr1);

llr2 =  abs(1/sqrt(2)+ img_rx_sig);   
llr3 =  abs(-1/sqrt(2)+ img_rx_sig);    
llr_quad = log(llr2./llr3);
% llr0 =  abs(1 + real_rx_sig);  
% llr1 =  abs(-1 + real_rx_sig);    
% llr_inphase = log(llr0./llr1);
% 
% llr2 =  abs(1+ img_rx_sig);   
% llr3 =  abs(-1+ img_rx_sig);    
% llr_quad = log(llr2./llr3);

n=0;
demod_output=zeros(1,G);
for k=1:1:length(Interleaved)/2
    n=n+1;
   demod_output(n)= llr_inphase(k); 
   n=n+1;
   demod_output(n)= llr_quad(k);
end
%--------------QPSK De-Modulation Done-----------%
%--------------De-Srambling--------------%

for d=1:length(demod_output)
    if x4(d)==0
        continue
    elseif x4(d)==1
        demod_output(d)=-demod_output(d);
    end
end
%--------------De-Srambling done--------------%

%-------------Segmentation--------%
X=reshape(demod_output, [E,C]);

%-----------------De-Interleaving----------------------%
Y=X;
Y=transpose(Y);
a1=E/order;
Z=zeros(a1,order);
Z1=zeros(C,E);
for h=1:C
    for h1=1:order
        Z(:,h1)=Y(h,h1:order:end);
    end
    Z1(h,:)=reshape(Z, 1, []);
end
Z1=transpose(Z1);
Z11=Z1;
%-----------------De-Interleaving Done----------------------%

%----------------------Rate Recovery-------------------%
for j=E+1:66*Zc
    Z11(j,:)=0.5;
end
%------------------Rate Recovery Done-----------------------%



%----------------LDPC decoding---------------------------%
outputbits= double(LDPCDecode(Z11,nbg,nldpcdecits));
% Every CB segment is LDPC Coded & Decoded and stored in 2-D matrix Outputbits  

outputbits=outputbits';


filler_removed=zeros(C,K1);
for i=1:K1
    filler_removed(:,i)=outputbits(:,i);
end
% FIller Bits have been removed
u=1;
%---------CB-CRC Validation Started---------%
for i=1:C
    if u==1
        remainder=crc_validate(filler_removed(i,:),generator);  
        for j=1:length(remainder)
            %for loop for checking if remainder is all 0's or
            %something else
            if remainder(j)==0
                ans1=1;
            else
                ans1=0;
                u=0;
                break;
            end
        end
    end
end
if u==1
    fprintf("%s\n","Code Block CRC validated Successfully");
else
    fprintf("%s\n","Code Block CRC validation Failed");
end

crc_removed=zeros(C,p);
%----CB-CRC Removed and Segmented Blocks stored in crc_removed matrix
for i=1:p
    crc_removed(:,i)=filler_removed(:,i);
end
crc_removed=transpose(crc_removed);

% Segmented Code Blocks are concatenated and stored in "final_output"
final_output = reshape(crc_removed, 1, []);

remainder=crc_validate(final_output,generator); %crc_validatation of TB-CRC
                  % is performed

for j=1:length(remainder) %for loop for checking if remainder is all 0's or
    %something else
    if remainder(j)==0
        ans2=1;
    else
        ans2=0;
        break
    end
end
if ans2==1
    fprintf("%s\n","Transport Block CRC validated Successfully");
else
    fprintf("%s\n","Transport Block CRC validation Failed");
end