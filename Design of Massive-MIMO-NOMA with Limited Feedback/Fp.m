function z = Fp(x,ai, N, Mt)

z = 1-gammainc( real(1/ai/x), N-Mt+1);
