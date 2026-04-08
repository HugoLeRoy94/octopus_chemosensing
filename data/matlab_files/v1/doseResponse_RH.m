function [hillCoeff, ec50, minDose, maxDose, coeffs, meanResponse, doses, sigmoid]=doseResponse_RH(dose,response, polarity)

% DOSERESPONSE     Computes the Hill Coefficient and EC50 of a
% dose/response relationship given two vectors describing the doses and
% responses. A semilog graph is plotted illustrating the relationship, upon
% which the man and standard error of the response to each dose level is
% plotted along with the fitted Hill Equation sigmoid. The EC50 is also
% labelled. Requires nlinfit from statistics toolbox.
%
%   Example
%       d=[3.75 3.75 3.75 3.75 3.75 7.5 7.5 7.5 7.5 15 15 15 15 15 60 60 60 60]
%       r=[107 91 99 124 100 96 92 133 119 84 66 86 106 91 52 37 10 69]
%       [hill ec50]=doseResponse(d,r)
%
%   Inputs
%       Two vectors of the same length, the first containing the dose and
%       the second the response. If doses of 0 are contained in the data,
%       these are used as control values and the data are normalised by
%       their mean value.
%
%   Notes
%       This function was written to rapidly produce simple dose/response curves
%       and EC50s for publication.

%deal with 0 dosage by using it to normalise the results.


tmpdose = dose;
tmpdose(tmpdose == 0) = 10^(log10(tmpdose(2))-2);

normalised=0;
if (sum(dose(:)==0)>0)
    %compute mean control response
    controlResponse=mean(response(dose==0));
    %remove controls from dose/response curve
    response1=response./controlResponse;
    response=response(dose~=0)./controlResponse;

    dose=dose(dose~=0);
    normalised=1;
end

%hill equation sigmoid
sigmoid=@(beta,x)beta(1)+(beta(2)-beta(1))./(1+(x/beta(3)).^beta(4));


%calculate some rough guesses for initial parameters
if polarity == -1
    minResponse=min(response1);
    maxResponse=max(response1);
    midResponse=mean([minResponse maxResponse]);
    minDose=min(tmpdose);
    maxDose=max(tmpdose);
elseif polarity == 1
    minResponse=max(response1);
    maxResponse=min(response1);
    midResponse=mean([minResponse maxResponse]);
    minDose=min(tmpdose);
    maxDose=max(tmpdose);
end


%fit the curve and compute the values
[coeffs,r,J]=nlinfit(tmpdose,response1,sigmoid,[minResponse maxResponse midResponse 1]);

ec50=coeffs(3);
hillCoeff=coeffs(4);

doses=unique(tmpdose);
meanResponse=zeros(1,length(doses));
stdErrResponse=zeros(1,length(doses));
for ss=1:length(doses)
    responses=response1(tmpdose==doses(ss));
    meanResponse(ss)=mean(responses);
    stdErrResponse(ss)=std(response1)/sqrt(length(response1));
    %stdResponse(ss)=std(responses);
end
