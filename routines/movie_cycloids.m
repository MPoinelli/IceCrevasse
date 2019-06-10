run tracking_rate.m

k = 1;

z = CRACK(k).TIME;
x = CRACK(k).XX;
y = CRACK(k).YY;
c = CRACK(k).PROPAGATION_RATES;

figure
curve = animatedline('LineWidth',2);
set(gca,'XLim',[min(x) max(x)],'YLim',[min(y) max(y)]);
hold on;grid on
for i=1:length(x)
    addpoints(curve,x(i),y(i)),z(i);
    head = scatter3(x(i),y(i),z(i),'filled','MarkerFaceColor','b');
    view(0,90)
    time = CRACK(k).TIME(i);
    string = {'T = ',num2str(time/(3600*24)),'days'};
    Text = text(-30,20,string);
    drawnow
    pause(0.01);
    delete(head);
    if i ~=length(x)
    delete (Text);
    end
end