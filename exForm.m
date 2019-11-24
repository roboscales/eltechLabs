%1.find a module of complex
%2.find a angel of complex in degrees
%3.print this in form: %4.3e^(j%4.3deg) fprintf()

function E  =  exForm(complex)
%display(E);

for c=1:length(complex)
    if angle(complex(c)) < 0
        E(c) = sprintf("%4.3f\\bullete^(-j%4.3f\\degree)", abs(complex(c)), abs(rad2deg(angle(complex(c)))));
    elseif angle(complex(c)) > 0
        E(c) = sprintf("%4.3f\\bullete^(j%4.3f\\degree)", abs(complex(c)), abs(rad2deg(angle(complex(c)))));     
    end
end
display('Eulers Form');
display(E');
end