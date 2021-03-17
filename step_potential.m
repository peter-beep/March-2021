function step_potential()
format longEng

time_vec = zeros(1,1000);

test_vec = linspace(1,500,500);
time_vec(1) = 0;

J=50;
matrix_test = zeros([J 500]);



for J = 1 : 4
    for AJ = 1 : 250
        for height = 0.004: 0.1 : 25.004
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
    end 
end 

    J=2;
    for AJ = 250 : 490
        for height = 0.003 : 0.05 : 0.015
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
    end
    
    for AJ = 501 : 700
        for height = 0.008 : -0.5 : 0.002
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
    end
    
    J=3;
     for AJ = 251 : 400
        for height = 0.007 : -0.5 : 0.001
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
     end
     
     disp('third row')
     
     for AJ = 601 : 900
        for height = 0.009:-0.5:0.002
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
     end
     

    J=4;
    for AJ = 251 : 800
        for height = 0.01:-0.5:0.002
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
    end
     
    for AJ = 801 : 900
        for height = 0.0009:-0.5:0.0004
            sigma = 0.4./height;
            matrix_test(J,AJ) = sigma;
        end      
     end

for next_temp = 1 : 46
    x = matrix_test(1,:);
    x = x(randperm(length(x)));
    matrix_test(4+next_temp,:)=x;
end 


temp_exit=0;
for B = 3: 50
    for A = 1 : length(matrix_test(B,:))
    % construct gaussians 
    mu = 0.99999999999+ (A-1);
    disp('mu')

    syms y(x);
    % initialize different potential from product of all distributions
    sigma_temp = matrix_test(B,A);
    height = 0.4./sigma_temp;
    func_exp_1 = @(x, mu, sigma) exp(-((x-(mu))./sigma).^2);
    func_next = func_exp_1(x , mu + A-1, height);
    
    % construct product of variances
    func_variance = @(sigma) sqrt(2 .* pi .* sigma);
    temp_next = func_variance(matrix_test(B,A));
    func_boltzmann_temp = @(temp) temp;
    % 5500, 2500, 1000
    bt_temp = func_boltzmann_temp(70000); 
    disp(bt_temp)
    
    % store result
    syms y(x);
    ode = diff(y,x,2) == - 1./(bt_temp) + (func_next) .* diff(y,x,1)./(bt_temp);
    Dy = diff(y);
    cond_1 = y(A)==0;
    cond_2 = Dy(0)==0;
    conds_temp = [cond_1 cond_2];
    ySol = dsolve(ode,conds_temp);

    pretty(ySol)
    y = sym(ySol);
    y = char(y);
    temp_sol=y;
    disp('sol')
    % process string contents of solution
    temp_pos = strfind(y,'int');
    ySol = dsolve(ode,conds_temp);
    
    pretty(ySol)
    y = sym(ySol);
    % temp_sol = y;
    % array_temp = zeros(1,7 * length(test_vec));
   
        temp_index=1;
            % array_temp(temp_index) = 
            disp('begin')
            disp('exponential term')
            temp_str_zero = temp_sol(1: temp_pos(temp_index)-2);
            disp('fist piece')
            disp(temp_str_zero)
            
            % first part of exponent 
            

            temp_str_one = temp_sol(temp_pos(temp_index) :temp_pos(temp_index+1)-3);
            disp('first integral term')
            disp(temp_str_one)
            x_1 = strfind(temp_str_one, '(');
            x_2 = strfind(temp_str_one,',');
            % temp_str_one = temp_sol(temp_pos(temp_index)+1 :temp_pos(temp_index+1)-3);
            % array_temp(temp_index+1) = 
            % temp_str_one = str2num(temp_sol(temp_pos(temp_index) :temp_pos(temp_index+1)-3));
            temp_str_one = temp_sol(temp_pos(temp_index) :temp_pos(temp_index+1)-3);
            % disp(temp_str_one)
            % loop over temp_str to numerically integrate part of solution
            x_1 = strfind(temp_str_one,'(');
            x_2 = strfind(temp_str_one,',');
            temp_str_one = temp_str_one(x_1(1)+1: x_2(1)-1);
            disp('string output')
            % disp(class(temp_str_one))
            disp(temp_str_one)
            disp(length(temp_str_one))
            % place all contents of string in first entry
            temp_str_one=convertCharsToStrings(temp_str_one);
            disp(class(temp_str_one))
            % disp(temp_str_one)
            % temp_str_one=str2double(temp_str_one);
            syms u;
            % temp_str_one = append("@(u)", temp_str_one);
            disp(temp_str_one)
            % temp_str_one = str2sym(char(temp_str_one));
            syms u;
            syms x;
            func_temp = matlabFunction(str2sym(temp_str_one));

            disp('FIRST, SECOND AND THIRD TERMS')
            
            trap_func_1 = sym(func_temp(0));
            disp(trap_func_1)
            trap_func_2 = sym(func_temp(u./2));
            disp(trap_func_2)
            trap_func_3 = sym(func_temp(u));
            disp(trap_func_3)
            trap_final = (u./2) .* symsum(trap_func_1 , trap_func_2,trap_func_3);
            trap_final = trap_final .* str2sym(temp_str_one);
            final_handle = convertCharsToStrings(trap_final);
            disp(trap_final)
            integral_term_four_five = vpaintegral(final_handle,A-1, A);

            final= @(x) -(A^2./bt_temp - x^2./bt_temp - integral_term_four_five);
            disp(final(A-1))
            tmpr = final(A-1);
            temp_exit = tmpr+temp_exit;
            
            disp('end')

            % exp((1329227995784915872903807060280344576*pi^(1/2)*erf(25000*u - 12500))/14546975136154065625) .* exp((1329227995784915872903807060280344576*pi^(1/2)*erf(25000*v - 12500))/14546975136154065625),  1, 2,0,1)


            figure(B+3)
            subplot(4,1,1)
            title('Relationship between exit time and base pair of the target sequence')
            xlabel('Base pair')
            ylabel('Exit time')
            plot(A,double(temp_exit),'.')
            hold on;
            
            % plot inverse of exit times
            subplot(4,1,2)
            title('Inverse of the exit time')
            xlabel('Base pair')
            ylabel('Reaction time')
            plot(A,-1./double(temp_exit),'.')
            hold on;
            
            subplot(4,1,3)
            title('Relationship between barrier height and magnitude of exit time')
            xlabel('Base pair')
            ylabel('Fluctuation height')
            plot(A,matrix_test(B,A),'.')
            hold on;
            
            subplot(4,1,4)
            title('Three dimensional solution curve')
            xlabel('Base pair')
            ylabel('Exit time')
            zlabel('Temperature')
            plot3(double(temp_exit), A, bt_temp,'.')
            hold on;
            
            bt_temp = bt_temp+10;
    end   
end  
end  
    
    % figure(temp);
    % plot(Cell, linspace(1,10,10))
    % hold on;
    
    % hold off;
 
 


