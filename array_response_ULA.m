function y = array_response_ULA(a1,N)
y = exp(-1i*pi*(0:N-1)'*sin(a1));
end