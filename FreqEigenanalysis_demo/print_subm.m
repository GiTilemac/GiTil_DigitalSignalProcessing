% Prints results of 'FreqEigenanalysis'
function [  ] = print_subm( subm, A, sigma_w, omega )
    N = numel(subm);
    if(N == 6)
        fprintf('Parameter M:               2              10              20              30              40              50\n');
    else
        fprintf('Parameter M:               10              20              30              40              50\n');
    end
    if(numel(omega)==1)
        fprintf('Amplitude: %.4f    ',A);
        for i = 1:N
            fprintf('%.4f(%.2f%%)  ',subm{i}.A,100*(1-abs(subm{i}.A-A)/A));
        end
    else
        for j = 1:numel(omega)
            fprintf('Amplitude: %.4f    ',A(j));
            for i = 1:N
                fprintf('%.4f(%.2f%%)  ',subm{i}.A(j),100*(1-min( abs(subm{i}.A(j)-A(j))/A(j) , abs(subm{i}.A(j)-A(j))/subm{i}.A(j))));
            end
            fprintf('\n');
        end 
    end
    fprintf('\n');
    fprintf('Noise variance: %.1f  ',sigma_w);
    for i = 1:N
        fprintf('%.2f(%.2f%%)    ',subm{i}.sigma,100*(1-min( abs(subm{i}.sigma-sigma_w)/sigma_w , abs(subm{i}.sigma-sigma_w)/subm{i}.sigma)));
    end
    fprintf('\n');
    if(numel(omega)==1)
        fprintf('Frequency: %.4f    ',omega);
        for i = 1:N
            fprintf('%.4f(%.2f%%)  ',subm{i}.omega,100*(1-abs(subm{i}.omega-omega)/omega));
        end
    else
        for j = 1:numel(omega)
            fprintf('Frequency: %.4f    ',omega(j));
            for i = 1:N
                fprintf('%.4f(%.2f%%)  ',subm{i}.omega(j),100*(1-min(abs(subm{i}.omega(j)-omega(j))/omega(j),abs(subm{i}.omega(j)-omega(j))/subm{i}.omega(j)) ));
            end
            fprintf('\n');
        end
    end
    fprintf('\n');
    fprintf('1st order central tendency: ');
    for i = 1:N
        fprintf('%.2e    ',subm{i}.cten1);
    end
    fprintf('\n');  
    fprintf('2nd order central tendency: ');
    for i = 1:N
        fprintf('%.2e    ',subm{i}.cten2);
    end
    fprintf('\n\n\n'); 
end

