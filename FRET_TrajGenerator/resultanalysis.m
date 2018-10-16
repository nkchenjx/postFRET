for x = 1:lxy
   for y = 1:lxy
    for repeat = 1:10
     k12(repeat) = (R{x,y,repeat}.rateconstantsfit4(2,1)/rateSeq{x,y}(2,1)-1)*100;
    end
        K12(x,y) = mean(k12);
    K12std(x,y) = std(k12);
end
end
figure; imagesc(K12); title('K12');
figure; imagesc(K12std); title('K12std');

for x = 1:lxy
   for y = 1:lxy
    for repeat = 1:10
     k12(repeat) = (R{x,y,repeat}.rateconstantsfit4(2,1)/R{x,y,repeat}.rateconstantsfit1(2,1)-1)*100;
    end
    K12t(x,y) = mean(k12);
    K12tstd(x,y) = std(k12);
end
end

figure; imagesc(K12t); title('K12t');
figure; imagesc(K12tstd); title('K12tstd');


for x = 1:lxy
   for y = 1:lxy
    for repeat = 1:10
     k21(repeat) = (R{x,y,repeat}.rateconstantsfit4(1,2)/rateSeq{x,y}(1,2)-1)*100;
    end
    K21(x,y) = mean(k21);
    K21std(x,y) = std(k21);
end
end
figure; imagesc(K21); title('K21');
figure; imagesc(K21std); title('K21std');

for x = 1:lxy
   for y = 1:lxy
    for repeat = 1:10
     k21(repeat) = (R{x,y,repeat}.rateconstantsfit4(1,2)/R{x,y,repeat}.rateconstantsfit1(1,2)-1)*100;
    end
    K21t(x,y) = mean(k21);
    K21tstd(x,y) = std(k21);
end
end

figure; imagesc(K21t); title('K21t');
figure; imagesc(K21tstd); title('K21tstd');



for x = 1:lxy
   for y = 1:lxy
    for repeat = 1:10
     k12(repeat) = (R{x,y,repeat}.rateconstantsfit2(2,1)/rateSeq{x,y}(2,1)-1)*100;
    end
    K12NF(x,y) = mean(k12);
    K12NFstd(x,y) = std(k12);
end
end
figure; imagesc(K12NF); title('K12NF');
figure; imagesc(K12NFstd); title('K12NFstd');
