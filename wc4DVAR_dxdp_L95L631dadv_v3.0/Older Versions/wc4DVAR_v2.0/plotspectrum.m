figure;
subplot(3,2,5)
plot(DSp), title('n=12'), xlabel('\lambda_i(S_p)'), ylabel('Value');
subplot(3,2,6)
plot(DSx,'r'), title('n=12'), xlabel('\lambda_i(S_x)'), ylabel('Value');

figure;
semilogy(DSp), title('n=12'), xlabel('\lambda_i(S_p)'), ylabel('Value');
hold;
semilogy(DPSp,'--'), 
figure;
semilogy(DSx,'r'), title('n=12'), xlabel('\lambda_i(S_x)'), ylabel('Value');
