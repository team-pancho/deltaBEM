function UG = unpackGeometry(MG)

% UG = unpackGeometry(MG)
% Input:
%       MG    :    merged geometry (one structure array)
% Output:
%       UG    :    cell array (several elements; one structure array per
%                  element)
%
% Last modified: May 14, 2014

K = length(MG.comp);
UG{K} = [];
MG.comp=[MG.comp length(MG.next)+1]; 

for i=1:K
    next=MG.next(MG.comp(i):MG.comp(i+1)-1)-MG.comp(i)+1; 
    UG{i} = struct('midpt',MG.midpt(MG.comp(i):MG.comp(i+1)-1,:),...
                   'brkpt',MG.brkpt(MG.comp(i):MG.comp(i+1)-1,:),...
                   'normal',MG.normal(MG.comp(i):MG.comp(i+1)-1,:),...
                   'next',next,...
                   'comp',1);
    if ~isfield(MG,'parity')
        continue
    elseif MG.parity(MG.comp(i),MG.comp(i))==1
            UG{i}.parity = speye(length(UG{i}.midpt));
            UG{i}.parity = UG{i}.parity-UG{i}.parity(end:-1:1,:);
    end 
end
return

