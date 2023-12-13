function [code_align] = Codetransform_HRL(align_seq,code)
%UNTITLED 转换code 从没有序列比对的原始序列上的code转换到经过序列比对后的code
%   align_seq是序列比对后包含gap(-)的序列
%   code是原始序列（没有gap）上的code位置
%   code_align是code在序列比对后的sequence上的位置
gap_site = strfind(align_seq,'-');
code_align = zeros(length(code),1);
for i = 1:length(code)
    new = code(i)+sum(gap_site<=code(i));
    old = code(i);
    while sum(gap_site<=old) ~= sum(gap_site<=new)
        old = new;
        new = code(i)+sum(gap_site<=new);
    end
    code_align(i) = new;
end
end

