clear()
x0 = -5;
xn = 5;
x = linspace(x0,xn,11);
x_org = linspace(x0,xn,101);
% y_org = sin(x_org);
y_org = zeros(1,101);
y_org(50) = 1;
y = zeros(1,11);
y(6) = 1;

A(:,1) = ones(length(x),1);
A(:,2) = x;
A(:,3) = x.^2;
A(:,4) = x.^3;
A(:,5) = x.^4;
A(:,6) = x.^5;
A(:,7) = x.^6;
A(:,8) = x.^7;
A(:,9) = x.^8;
A(:,10) = x.^9;

b_hat = inv(A'*A)*A'*y';


    

figure
plot(x,y,'o', 'DisplayName', 'Fit Points')
hold on
plot(x, y, 'b','LineWidth', 2,'DisplayName', 'og func.');
hold on
    
for i= 1
    
    pol = @(x) b_hat(1).*ones(1,length(x)) + b_hat(2).*x + b_hat(3).*x.^2 ...
    +b_hat(4).*x.^3+b_hat(5).*x.^4 +b_hat(6).*x.^5+b_hat(7).*x.^6 ...
    +b_hat(8).*x.^7 %+ b_hat(9).*x.^8 + b_hat(10).*x.^9;

    y_hat = pol(x_org);
    plot(x_org,y_hat, '--','DisplayName',['N = ' num2str(10-i)] );
    axis([x0 xn -0.5 1.5]);
    b_hat(length(b_hat)-i+1) = 0;
    title("Non-Symmetric Peak");
%     b_hat

end
legend(gca,'show')

