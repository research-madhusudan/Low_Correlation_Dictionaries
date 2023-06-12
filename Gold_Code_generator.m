close all
clear all
clc
tic
% Register length : L=5,6,7,9,10,11
L = 10;
% Feedback taps (generator polynomial) in descending order
if mod(L,4)==0
    display("ERROR: gold-code doesn't exist when L is a multiple of 4")
    display("Try again with different value of L. Killing the Program.....")
    return
elseif L==5
    feedbackTaps1=[5 2];
    feedbackTaps2=[5 4 3 2];
elseif L==6
    feedbackTaps1=[6 1];
    feedbackTaps2=[6 5 2 1];
elseif L==7
    feedbackTaps1 = [7 3];
    feedbackTaps2 = [7 3 2 1];
elseif L==9
    feedbackTaps1=[9 4];
    feedbackTaps2=[9 6 4 3];
elseif L==10
    feedbackTaps1=[10 3];
    feedbackTaps2=[10 8 3 2];
elseif L==11
    feedbackTaps1=[11 2];
    feedbackTaps2=[11 8 5 2];
end

% length of m-sequences
lenSequence = 2^L - 1;

% Store the generated m-sequences
mSequence1 = zeros(lenSequence,1);
mSequence2 = zeros(lenSequence,1);

% initial condition to generate m-sequence
initialCondition = int2bit(1, L)';

% Initialize the shift register with the current initial condition
shiftRegister1 = initialCondition;
shiftRegister2 = initialCondition;

% Generate the m-sequence
for i = 1:lenSequence
    mSequence1(i, 1) = shiftRegister1(end);
    mSequence2(i, 1) = shiftRegister2(end);

    feedback1 = mod(sum(shiftRegister1(feedbackTaps1)), 2);
    feedback2 = mod(sum(shiftRegister2(feedbackTaps2)), 2);

    shiftRegister1 = [feedback1 shiftRegister1(1:L-1)];
    shiftRegister2 = [feedback2 shiftRegister2(1:L-1)];
end

% initialize set of Gold sequences using original m-sequences 
GoldSequences=[mSequence1,mSequence2];

% populate the set of Gold sequences using remaining sequences of the form XOR(mSequence1,(T^i)*mSequence2)
I=eye(lenSequence);
for i=0:lenSequence-1
    GoldSequences=[GoldSequences mod(mSequence1+circshift(I,i)*mSequence2,2)];
end
B=[];
for i=1:lenSequence
        B=[B circshift(I,i)*GoldSequences];
end

% convert {0,1} sequence to {-1,+1} sequence and normalize the columns to
% unit norm
B=2*B-1;
% uncomment following lines to check the correlation properties of the
% dictionary matrix formed using all shift values of Gold Codes
% unique(B'*B)'
% imagesc(abs(B'*B))
B=normc(B);

% add identity matrix to obtain larger dictionay matric with same mutual
% coherence properties
I=eye(lenSequence);
B=[B I];
txt1="Gold_sequences_size_"+string(lenSequence)+".mat";
txt2="Gold_based_dictionary_"+string(lenSequence)+"x"+string(size(B,2))+".mat"
save(txt1,"GoldSequences")
save(txt2,"B")
toc