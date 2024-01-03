function  setSMTP()
%SETSMTP This function set the email smtp sender
EMAIL_AUTHCODE = "SBYEKBGMPXXNUELU";
EMAIL_SEND = "liweihanqq@163.com";

setpref('Internet', 'E_mail', EMAIL_SEND);
setpref('Internet', 'SMTP_Server', 'smtp.163.com');
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth', 'true');
setpref('Internet', 'SMTP_Username', EMAIL_SEND);
setpref('Internet', 'SMTP_Password', EMAIL_AUTHCODE);
end