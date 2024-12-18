function ovsf_codes1=ovsf_generator(spread_factor,code_number)
ovsf_codes=1;
for i=1:1:log2(spread_factor)
      temp=ovsf_codes;
  for j=1:1:size(ovsf_codes,1)
     if j==1
          ovsf_codes=[temp(j,:),temp(j,:); temp(j,:),(-1)*temp(j,:)];
     else
          ovsf_codes=[ovsf_codes; temp(j,:),temp(j,:); temp(j,:),(-1)*temp(j,:)];
    end
  end
end
ovsf_codes1=ovsf_codes(code_number,:);%OVSFÂëÉú³É