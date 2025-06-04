function obj = aesobj()
%AESOBJ This function use a fix code as key and generate the aes object
% The file will be generate as pcode for hidding the details

% public key: "REGISTER3D_KEY"
% algorithm: SHA-1

obj = AES("REGISTER3D_KEY", "SHA-1");
end