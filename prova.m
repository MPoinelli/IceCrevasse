pizzi

NO = [100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000]

for oba =1:length(NO)
    Nodes = NO(oba)
    run smooth_Europa_features.m
    run tracking_rate.m
    
    abc(oba) = mean(ACTIVITY);
    def(oba) = median(ACTIVITY);
    clear ACTIVITY
end

figure,semilogx(NO,abc)
figure,semilogx(NO,def)