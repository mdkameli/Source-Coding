clc; clear all; close all;

%% Read Input Audio signal

%%% Classic Music %%%
%Audio = 'chopin.wav';

%%% Jazz Music %%%
%Audio = 'Jazz-Test.wav';

%%% Traditional Voice %%%
%Audio = 'Traditional.wav';

%%% POP Music %%%
%Audio = 'POP.wav';

%%% Folk Music %%%
%Audio = 'Folk.wav';

%%% Woman Voice %%%
Audio = 'Woman.wav';

%%% Man Voice %%%
%Audio = 'Man.wav';

[x,Fs] = audioread(Audio,'native');                          %Read Input Audio
l = length (x);
fr = (0:l-1)*(Fs/l);

%% Calculate Prediction Sequence Parameter
%%% Calculation Aoutocorrelation , R
%%% Order N =1 

dbl = double(x); 
v = sum(dbl.^2)/l;
R(1) = sum(dbl(1:l-1).*dbl(2:l))/(l-1);
R(2) = sum(dbl(1:l-2).*dbl(3:l))/(l-2);
R(3) = sum(dbl(1:l-3).*dbl(4:l))/(l-3);
R(4) = sum(dbl(1:l-4).*dbl(5:l))/(l-4);


%%% Calculation of Pn Coefficient
%%% For Orders of N = 1, 2 & 4
%%% a1 coefficient order 1
%%% a2 coefficient order 2
%%% a4 coefficient order 4

%%% Order N = 1
a1 = R(1)/v;                                %Coefficient Order 1

%%% Order N = 2
Z2 = [v R(1); R(1) v];
R2 = [R(1);R(2)];                           %Coefficient Order 2
a2 = Z2\R2;

%%% Order N =4
Z4 = [v R(1) R(2) R(3);R(1) v R(1) R(2);R(2) R(1) v R(1);R(3) R(2) R(1) v];
a4 = Z4\R';                                 %Coefficient Order 4

%% DPCM (Differential Pulse-Code Modulation Implementation
%%% Calculate Differential Argument and Quantize data

M = 64;                                           % Number of bit 
%M = 128;                                           % Number of bit 
%M = 256;                                           % Number of bit 
%M = 2048;
B = (2^16);                                        % Total Width  
delta = (B/M);                                     % Decision width

%%% Order N =1
%x = double(x);
d1(1) = x(1);
xhat1(1) = x(1); 
 for k = 2:l
    d1(k) = x(k) - a1*xhat1(k-1);
    if d1(k) >= (B/2) 
        dhat1(k) = (B/2) + (delta/2);
    elseif d1(k) <= -(B/2)
        dhat1(k) = -(B/2)- (delta/2);
    else 
        dhat1(k) = delta*floor((d1(k)/delta));    
    end
    xhat1(k) = a1*xhat1(k-1) + dhat1(k);
end
% filename = 'Man-Order1-64Level.wav';
% audiowrite(filename,xhat1,Fs);
vq1 = sum((x-xhat1').^2)/l;
SNR1 = 10*log(v/vq1);
distor1 = 10*log(sum((x-xhat1').^2)/l);        % Mean square error

%%% Order N =2

d2(1) = x(1);
xhat2(1:2) = x(1:2);                           %Initial Values fo Predictor
for k = 3:l
    d2(k) = x(k) - a2(1)*xhat2(k-1) - a2(2)*xhat2(k-2);
    if d2(k) >= (B/2) 
        dhat2(k) = (B/2)+ delta/2;
    elseif d2(k) <= -(B/2)
        dhat2(k) = -(B/2)- delta/2;
    else 
        dhat2(k) = delta*floor((d2(k)/delta));
    end
    xhat2(k) = a2(1)*xhat2(k-1) + a2(2)*xhat2(k-2) + dhat2(k);
end
% filename = 'Man-Order2-64Level.wav';
% audiowrite(filename,xhat2,Fs);                      %Write Output File
vq2 = sum((x-xhat2').^2)/l;
SNR2 = 10*log(v/vq2);                                 %Calculate SNR
distor2 = 10*log(sum((x-xhat2').^2)/l);               %Calculate Distortion
 
%%% Order N =4

d4(1) = x(1);                                  %Initial Values fo Predictor
xhat4(1:4) = x(1:4);                                
for k = 5:l
    d4(k) = x(k) -a4(1)*xhat4(k-1) -a4(2)*xhat4(k-2) -a4(3)*xhat4(k-3) -a4(4)*xhat4(k-4);
    if d4(k) >= (B/2) 
        dhat4(k) = (B/2)+ delta/2;
    elseif d4(k) <= -(B/2)
        dhat4(k) = -(B/2)- delta/2;
    else 
        dhat4(k) = delta*floor((d4(k)/delta));
    end
    xhat4(k) = a4(1)*xhat4(k-1) +a4(2)*xhat4(k-2) +a4(3)*xhat4(k-3) +a4(4)*xhat4(k-4) + dhat4(k);
end
% filename = 'Man-Order4-64Level.wav';
% audiowrite(filename,xhat4,Fs);                      %Write Output File
vq4 = sum((x-xhat4').^2)/l;
SNR4 = 10*log(v/vq4);               %Calculate SNR
distor4 = 10*log(sum((x-xhat4').^2)/l);                     %Calculate Distortion

%% Mapping Method
%%% Map the Prediction error to positive Integer and Normalize
e = [];
for k = 1:100000
      if dhat1(k) < 0
        e = [e -dhat1(k)-(delta/2)];
      else
        e = [e dhat1(k)];
      end
end

%%% Geometric Distribution Parameter
clear p;clear m;
p1 = length(e)/(sum(e)+length(e));
m1 = floor(-1/log2(1-p1));

%%% Normalize Data
en = e./(delta/2);
%%% m Parameter
%%% First Approach
p2 = nnz(en)/length(en);
m2 = ceil(-1/log2(p2));
%%% Second approach
p3 = length(en)/(sum(en)+length(en));
m3 = ceil(-1/log2(1-p3));
%%
% indexa = find(dhat1>=0);
% indexb = find(dhat1<0);
% Pos = dhat1(dhat1>=0);
% Neg_mapped = dhat1(dhat1<0);
% Neg=-(Neg_mapped+(delta/2));
% En = [];
% for i=1:length(Pos)
%     En(Pos(i)) = indexa(i);
% end
% for i=1:length(Neg)
%     En(Neg(i)) = indexb(i);
% end

%% Golomb Coding

clear code;
clear q;
clear r;
code = [];
for i = 1:length(en);
    
    m = 2;
    n = en(i);
    q = floor(n/m);                   %compute the integer part   
    r = rem(n,m);                     %compute the reamind Part
    q1 = ones(1,q);       
    q_code = [q1 0];                  %unary code of quotient q
    
    [f,e]=log2(m);                    %f,e used to check if m is a power of 2
    if f==0.5 && e == 1               %special case of m=1
        code = q_code;
    else if f==0.5                    %check whether m is a power of 2                        
            r_code = de2bi(r,log2(m),'left-msb'); %log2(m)-bit binary code of r
        else
            a = ceil(log2(m));
            b = floor(log2(m));
        
            if r < (2^a - m)                
                r_code = de2bi(r,b,'left-msb'); %b-bit binary representation of r
            else
                r = r+(2^a - m);
                r_code = de2bi(r,a,'left-msb'); %a-bit binary representation of r
            end
        end
    end
newcode =[q_code r_code];               %golomb code of the input
code = [code newcode];                  %Binary Sequence of Codes
end
R1 = length(code)/length(en);

%% Decoder
%%% Golomb Decoder
code_len = length(newcode);
%%% Count the number of 1's followed by the first 0
q = 0;
for i = 1: code_len
    if newcode(i) == 1
        q = q + 1; %count the number of 1's
    else
        ptr = i;   % first 0
        break;
    end
end

if (m == 1)
    n_dec = q;   % special case for m = 1
else
    A = ceil(log2(m));
    B = floor(log2(m));
%%% decoding the remainder
    bcode = newcode((ptr+1): (ptr + B));
    r = bi2de(bcode,'left-msb');
    if r < (2^A - m)
        ptr = ptr + B;
    else
%%% r is the A-bit represtation of (r + (2^A - m))
        bcode = newcode((ptr+1): (ptr + A));
        r = bi2de(bcode,'left-msb') - (2^A - m);
        ptr = ptr + A;
    end
    n_dec = q * m + r; %computing the symbol from the decoded quotient and remainder
end

if ~isequal(ptr, code_len)
    error('Error: More than one codeword detected!');
end

%%% Map The Coded Data ("en" to "er")
%%% Denormalization and extracting dhat

for k = 1:length (en)
    if mod(en(k),2) == 1
        er(k) = -((en(k)+1)*(delta/2));
    else
        er(k) = en(k)*(delta/2);
    end
end

%%% Extracting Output
y1(1) = x(1);
for i = 2:length(er)
    y1(i) = a1*y1(i-1) + er(i);
end
 
%% Plot The Original Signal & Extracted Data
plot(fr,x,fr,xhat1,'--')
legend('Original signal','Decoded signal','Location','NorthOutside');

%%
%filename = 'Audio-Order1-8Level.wav';
%audiowrite(filename,xhat1,Fs);