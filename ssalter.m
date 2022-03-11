% s is secondary structure from Stride

s = '  TTTT  EEEEEETTEEEEE       EEEEEEEEEE       EEEEETTTTTEEEEE  TTTTEEEEEEEEEETTTTT      EEEEEE ';
prot = pdbread('fib3.pdb');
[res,n] = unique([prot.Model.Atom.resSeq]);
chain = [prot.Model.Atom.chainID];
chain = chain(n);
% alter /,ss='L'

for i = 1:length(s)
    if s(i)=='H'
        fprintf('alter chain %s and %g/,ss=''H''\n',chain(i),res(i))
    end
    if s(i)=='E'
        fprintf('alter chain %s and %g/,ss=''S''\n',chain(i),res(i))
    end
end

%rebuild