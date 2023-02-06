#include <cstdio>
#include <string>
#include <complex>
#include <cmath>
#include <unordered_map>
#include <set>
#include <memory>
#include <cassert>
#include <vector>
#include <functional>
#include <fftw3.h>
#include <list>
#include <unordered_set>

#include "fftgen.hpp"

int fft_nops_real[MAXFFT + 1];

int fft_nops[MAXFFT + 1];

void print_fft(int N) {
	std::string fname = "fft." + std::to_string(N) + ".cpp";
	set_file(fname);
	print("#include \"fft.hpp\"\n");
	print("\nvoid fft_%i(std::complex<double>* x0) {\n", N);
	indent();
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	print("double* x = reinterpret_cast<double*>(x0);\n");
	fft_bitreverse(N);
	fft(N, 0, true);
	fft_nops[N] = fft_opcnt(N, true);
	deindent();
	print("}\n\n");
}

void print_fft_real(int N) {
	std::string fname = "fft.real." + std::to_string(N) + ".cpp";
	set_file(fname);
	print("#include \"fft.hpp\"\n");
	print("\nvoid fft_real_%i(double* x0, double* y) {\n", N);
	indent();
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	print("double* x = reinterpret_cast<double*>(x0);\n");
	fft_bitreverse_real(N);
	fft_real(N, 0);
	print("y[%i] = x[%i];\n", 0, 0);
	print("y[%i] = 0;\n", 1);
	for (int n = 1; n < (N + 1) / 2; n++) {
		print("y[%i] = x[%i];\n", 2 * n, n);
		print("y[%i] = x[%i];\n", 2 * n + 1, N - n);
	}
	if (N % 2 == 0) {
		print("y[%i] = x[%i];\n", 2 * (N / 2), N / 2);
		print("y[%i] = 0;\n", 2 * (N / 2) + 1);
	}
	fft_nops_real[N] = fft_real_opcnt(N, true);
	deindent();
	print("}\n\n");
}

void print_fft_real_inv(int N) {
	std::string fname = "fft.real.inv." + std::to_string(N) + ".cpp";
	set_file(fname);
	print("#include \"fft.hpp\"\n");
	print("\nvoid fft_real_inv_%i(double* x0, std::complex<double>* y) {\n", N);
	indent();
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	print("double* x = reinterpret_cast<double*>(x0);\n");
	print("x[%i] = y[%i].real();\n", 0, 0);
	for (int n = 1; n < (N + 1) / 2; n++) {
		print("x[%i] = y[%i].real();\n", n, n);
		print("x[%i] = y[%i].imag();\n", N - n, n);
	}
	if (N % 2 == 0) {
		print("x[%i] = y[%i].real();\n", N / 2, N / 2);
	}
	fft_real_inv(N, 0);
	const auto J = fft_bitreverse_indices_real_inv(N);
	std::vector<int> I(N);
	for( int n = 0; n < N; n++) {
		I[J[n]] = n;
	}
	fft_bitreverse_real(N, I);
//	fft_nops_real[N] = fft_real_opcnt(N, true);
	deindent();
	print("}\n\n");
}

int main(int argc, char **argv) {

	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print_fft(n);
	}

	for (int n = 2; n <= MAXFFT; n++) {
		print_fft_real(n);
	}

	for (int n = 2; n <= MAXFFT; n++) {
		print_fft_real_inv(n);
	}

	set_file("fft.hpp");
	print("#include <complex>\n");
	print("#include <cmath>\n");
	print("#include <chrono>\n");
	print("#include <unordered_map>\n");
	print("#include <vector>\n\n");
	print("#include <cstring>\n\n");
	print("#include <fftw3.h>\n\n");
	print("#include <fenv.h>\n\n");
	print("\nvoid FFT(std::complex<double>* x, int N);\n");
	print("\nvoid FFT_real(std::complex<double>* xout, double* xin, int N);\n");
	print("\nvoid FFT_inv(std::complex<double>* x, int N);\n");

	std::string header4 = "inline std::complex<double> Ix(const std::complex<double>& a) {\n"
			"\tstd::complex<double> b;\n"
			"\tb.real(-a.imag());\n"
			"\tb.imag(a.real());\n"
			"\treturn b;\n"
			"}\n"
			"\n"
			"inline bool is_prime(int n) {\n"
			"\tbool v = true;\n"
			"\tif (n == 1) {\n"
			"\t\tv = false;\n"
			"\t} else {\n"
			"\t\tint kmax = sqrt(n);\n"
			"\t\tfor (auto i = 2; i <= kmax; i++) {\n"
			"\t\t\tif (n % i == 0) {\n"
			"\t\t\t\tv = false;\n"
			"\t\t\t\tbreak;\n"
			"\t\t\t}\n"
			"\t\t\tif (i > kmax) {\n"
			"\t\t\t\tbreak;\n"
			"\t\t\t}\n"
			"\t\t}\n"
			"\t}\n"
			"\treturn v;\n"
			"}\n"
			"\n"
			"inline int greatest_prime_factor(int N) {\n"
			"\tint v;\n"
			"\tfor (int n = 2; n <= N; n++) {\n"
			"\t\tif (N % n == 0 && is_prime(n)) {\n"
			"\t\t\tv = n;\n"
			"\t\t}\n"
			"\t}\n"
			"\treturn v;\n"
			"}"
			"\n"
			"inline bool allow_real(int N) {"
			"\treturn true;"
			"}";

	std::string header5 = "inline std::complex<double> No8x(const std::complex<double>& a) {\n"
			"\tstd::complex<double> c;\n"
			"\tstatic constexpr double c0 = sqrt(2) * 0.5;\n"
			"\tc.real(c0 * (a.real() + a.imag()));\n"
			"\tc.imag(c0 * (a.imag() - a.real()));\n"
			"\treturn c;\n"
			"}\n"
			"";
	std::string header41 = "inline std::complex<double> nIx(const std::complex<double>& a) {\n"
			"\tstd::complex<double> b;\n"
			"\tb.real(a.imag());\n"
			"\tb.imag(-a.real());\n"
			"\treturn b;\n"
			"}\n"
			"";

	std::string header51 = "inline std::complex<double> nNo8x(const std::complex<double>& a) {\n"
			"\tstd::complex<double> c;\n"
			"\tstatic constexpr double c0 = -sqrt(2) * 0.5;\n"
			"\tc.real(c0 * (a.real() + a.imag()));\n"
			"\tc.imag(c0 * (a.imag() - a.real()));\n"
			"\treturn c;\n"
			"}\n"
			"";

	std::string header42 = "inline std::complex<double> No8xIx(const std::complex<double>& a) {\n"
			"\tstd::complex<double> c;\n"
			"\tstatic constexpr double c0 = sqrt(2) * 0.5;\n"
			"\tc.imag(c0 * (a.real() + a.imag()));\n"
			"\tc.real(c0 * (a.real() - a.imag()));\n"
			"\treturn c;\n"
			"}\n"
			"";

	std::string header52 = "inline std::complex<double> nNo8xIx(const std::complex<double>& a) {\n"
			"\tstd::complex<double> c;\n"
			"\tstatic constexpr double c0 = -sqrt(2) * 0.5;\n"
			"\tc.imag(c0 * (a.real() + a.imag()));\n"
			"\tc.real(c0 * (a.real() - a.imag()));\n"
			"\treturn c;\n"
			"}\n"
			"";

	print("%s\n", header4.c_str());
	print("%s\n", header5.c_str());
	print("%s\n", header41.c_str());
	print("%s\n", header51.c_str());
	print("%s\n", header42.c_str());
	print("%s\n", header52.c_str());

	set_file("fft.cpp");
	print("#include \"fft.hpp\"\n\n");

	print("\n"
			"#include <chrono>\n"
			"\n"
			"class timer {\n"
			"\tstd::chrono::time_point<std::chrono::high_resolution_clock> start_time;\n"
			"\tdouble time;\n"
			"public:\n"
			"\tinline timer() {\n"
			"\t\ttime = 0.0;\n"
			"\t}\n"
			"\tinline void stop() {\n"
			"\t\tstd::chrono::time_point<std::chrono::high_resolution_clock> stop_time = std::chrono::high_resolution_clock::now();\n"
			"\t\tstd::chrono::duration<double> dur = stop_time - start_time;\n"
			"\t\ttime += dur.count();\n"
			"\t}\n"
			"\tinline void start() {\n"
			"\t\tstart_time = std::chrono::high_resolution_clock::now();\n"
			"\t}\n"
			"\tinline void reset() {\n"
			"\t\ttime = 0.0;\n"
			"\t}\n"
			"\tinline double read() {\n"
			"\t\treturn time;\n"
			"\t}\n"
			"};\n"
			"\n"
			"");

	print("\n");
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print("void fft_%i(std::complex<double>*);\n", n);
	}
	for (int n = 2; n <= MAXFFT; n++) {
		print("void fft_real_%i(double*, double*);\n", n);
	}
	for (int n = 2; n <= MAXFFT; n++) {
		print("void fft_real_inv_%i(double*, std::complex<double>*);\n", n);
	}
	print("\n");

	print("typedef void (*func_real_type)(double*, double*);\n");
	print("typedef void (*func_real_inv_type)(double*, std::complex<double>*);\n");
	print("typedef void (*func_type)(std::complex<double>*);\n");
	print("\n");
	print("const func_type fptr[] = {nullptr, nullptr, ");
	for (int n = 2; n <= MAXFFT; n++) {
		print_notab("&fft_%i", n);
		if (n != MAXFFT) {
			print_notab(", ");
		}
	}
	print_notab("};\n\n");
	print("const func_real_type fptr_real[] = {nullptr, nullptr, ");
	for (int n = 2; n <= MAXFFT; n++) {
		print_notab("&fft_real_%i", n);
		if (n != MAXFFT) {
			print_notab(", ");
		}
	}
	print_notab("};\n\n");
	print("const func_real_inv_type fptr_real_inv[] = {nullptr, nullptr, ");
	for (int n = 2; n <= MAXFFT; n++) {
		print_notab("&fft_real_inv_%i", n);
		if (n != MAXFFT) {
			print_notab(", ");
		}
	}
	print_notab("};\n\n");
	print("const int fft_nops[] = {0, 0, ");
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print_notab("%i", fft_nops[n]);
		if (n != MAXFFT) {
			print_notab(", ");
		}
	}
	print_notab("};\n\n");
	print("const int fft_nops_real[] = {0, 0, ");
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print_notab("%i", fft_nops_real[n]);
		if (n != MAXFFT) {
			print_notab(", ");
		}
	}
	print_notab("};\n\n");
	print("void FFT_real_inv(std::complex<double>* xout, double* xin, int N) {\n");
	indent();
	print("(*(fptr_real_inv[N]))(xin, xout);\n");
	deindent();
	print("}\n\n");
	print("void FFT_real(std::complex<double>* xout, double* xin, int N) {\n");
	indent();
	print("(*(fptr_real[N]))(xin, reinterpret_cast<double*>(xout));\n");
	deindent();
	print("}\n\n");
	print("void FFT(std::complex<double>* x, int N) {\n");
	indent();
	print("(*(fptr[N]))(x);\n");
	deindent();
	print("}\n\n");
	print("void FFT_inv(std::complex<double>* x, int N) {\n");
	indent();
	print("const double Ninv = 1.0 / N;\n");
	print("for( int n = 0; n < N; n++) {\n");
	indent();
	print("x[n] = std::conj(x[n]);\n");
	deindent();
	print("}\n");
	print("(*(fptr[N]))(x);\n");
	print("for( int n = 0; n < N; n++) {\n");
	indent();
	print("x[n] = std::conj(x[n]) * Ninv;\n");
	deindent();
	print("}\n");
	deindent();
	print("}\n\n");

	std::string header2 = "\n"
			"double rand1() {\n"
			"\treturn (rand() + 0.5) / RAND_MAX;\n"
			"}\n"
			"\n"
			"";
	std::string header3 = "\n"
			"void fftw(std::vector<std::complex<double>>& x) {\n"
			"\tconst int N = x.size();\n"
			"\tstatic std::unordered_map<int, fftw_plan> plans;\n"
			"\tstatic std::unordered_map<int, fftw_complex*> in;\n"
			"\tstatic std::unordered_map<int, fftw_complex*> out;\n"
			"\tif (plans.find(N) == plans.end()) {\n"
			"\t\tin[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);\n"
			"\t\tout[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * N);\n"
			"\t\tplans[N] = fftw_plan_dft_1d(N, in[N], out[N], FFTW_FORWARD, FFTW_ESTIMATE);\n"
			"\t}\n"
			"\tauto* i = in[N];\n"
			"\tauto* o = out[N];\n"
			"\tfor (int n = 0; n < N; n++) {\n"
			"\t\ti[n][0] = x[n].real();\n"
			"\t\ti[n][1] = x[n].imag();\n"
			"\t}\n"
			"\tfftw_execute(plans[N]);\n"
			"\tfor (int n = 0; n < N; n++) {\n"
			"\t\tx[n].real(o[n][0]);\n"
			"\t\tx[n].imag(o[n][1]);\n"
			"\t}\n"
			"\n"
			"}\n"
			""
			"void fftw_real(std::vector<std::complex<double>>& xout, const std::vector<double>& xin) {\n"
			"\tconst int N = xin.size();\n"
			"\tstatic std::unordered_map<int, fftw_plan> plans;\n"
			"\tstatic std::unordered_map<int, double*> in;\n"
			"\tstatic std::unordered_map<int, fftw_complex*> out;\n"
			"\tif (plans.find(N) == plans.end()) {\n"
			"\t\tin[N] = (double*) malloc(sizeof(double) * N);\n"
			"\t\tout[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));\n"
			"\t\tplans[N] = fftw_plan_dft_r2c_1d(N, in[N], out[N], FFTW_ESTIMATE);\n"
			"\t}\n"
			"\tauto* i = in[N];\n"
			"\tauto* o = out[N];\n"
			"\tfor (int n = 0; n < N; n++) {\n"
			"\t\ti[n] = xin[n];\n"
			"\t}\n"
			"\tfftw_execute(plans[N]);\n"
			"\tfor (int n = 0; n < N / 2 + 1; n++) {\n"
			"\t\txout[n].real(o[n][0]);\n"
			"\t\txout[n].imag(o[n][1]);\n"
			"\t}\n"
			"\n"
			"}\n"
			""
			"void fftw_real_inv(const std::vector<std::complex<double>>& xout, std::vector<double>& xin) {\n"
			"\tconst int N = xin.size();\n"
			"\tstatic std::unordered_map<int, fftw_plan> plans;\n"
			"\tstatic std::unordered_map<int, double*> in;\n"
			"\tstatic std::unordered_map<int, fftw_complex*> out;\n"
			"\tif (plans.find(N) == plans.end()) {\n"
			"\t\tin[N] = (double*) malloc(sizeof(double) * N);\n"
			"\t\tout[N] = (fftw_complex*) malloc(sizeof(fftw_complex) * (N / 2 + 1));\n"
			"\t\tplans[N] = fftw_plan_dft_c2r_1d(N, out[N], in[N], FFTW_ESTIMATE);\n"
			"\t}\n"
			"\tauto* i = in[N];\n"
			"\tauto* o = out[N];\n"
			"\tfor (int n = 0; n < N / 2 + 1; n++) {\n"
			"\t\to[n][0] = xout[n].real();\n"
			"\t\to[n][1] = xout[n].imag();\n"
			"\t}\n"
			"\t\to[0][1] = 0;\n"
			"\t\tif( N % 2 == 0 ) {\n"
			"\t\t\to[N / 2][1] = 0;\n"
			"\t\t}\n"
			"\tfftw_execute(plans[N]);\n"
			"\tfor (int n = 0; n < N; n++) {\n"
			"\t\txin[n] = i[n];\n"
			"\t}\n"
			"\n"
			"}\n"
			"";
	print("%s\n", header3.c_str());
	print("%s\n", header2.c_str());

	print("void test() {\n");
	indent();
	print("timer tm3, tm4;\n");
	print("for( int N = 2; N <= %i; N+=%i) {\n", MAXFFT, DFFT);
	indent();

	print_notab("\t\ttimer tm1, tm2;\n"
			"\t\tdouble err;\n"
			"\t\tdouble max;\n"
			"\t\tfor (int ti = 0; ti < 256; ti++) {\n"
			"\t\t\terr = 0.0;\n"
			"\t\t\tmax = 0.0;\n"
			"\t\t\tstd::vector<std::complex<double>> X(N);\n"
			"\t\t\tstd::vector<std::complex<double>> Y(N);\n"
			"\t\t\tfor (int n = 0; n < N; n++) {\n"
			"\t\t\t\tX[n] = std::complex<double>(rand1(), rand1());\n"
			"\t\t\t}\n"
			"\t\t\tauto X0 = X;\n"
			"\t\t\tY = X;\n"
			"\t\t\ttm1.start();\n"
			"\t\t\ttm3.start();\n"
			"\t\t\tFFT(X.data(), N);\n"
			"\t\t\ttm1.stop();\n"
			"\t\t\ttm3.stop();\n"
			"\t\t\ttm2.start();\n"
			"\t\t\ttm4.start();\n"
			"\t\t\tfftw(Y);\n"
			"\t\t\ttm2.stop();\n"
			"\t\t\ttm4.stop();\n"
			"\t\t\tfor (int i = 0; i < X.size(); i++) {\n"
			"\t\t\t\tY[i] -= X[i];\n"
			"\t\t\t}\n"
			"\t\t\tfor (int n = 0; n < N; n++) {\n"
			"\t\t\t\terr += std::abs(Y[n]) * std::abs(Y[n]);\n"
			"\t\t\t\tmax = std::max(max, std::abs(X0[n]));\n"
			"\t\t\t\t//printf(\"%i %e %e %e %e\\n\", n, X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag());\n"
			"\t\t\t}\n"
			"\t\t\terr = sqrt(err / N) / max;\n"
			"\t\t}\n"
			"\t\tprintf(\"%4i %4i %e %e %e %e %e %e %e\\n\", N, fft_nops[N], err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());\n"
			"");
	deindent();
	print("}\n/*return ;*/\n");
	print("for( int N = 2; N <= %i; N+=%i) {\n",  MAXFFT, DFFT);
	indent();
	print("if( !allow_real(N) ) {\n");
	indent();
	print("continue;\n");
	deindent();
	print("}\n");
	print_notab("\t\ttimer tm1, tm2;\n"
			"\t\tdouble err;\n"
			"\t\tdouble max;\n"
			"\t\tfor (int ti = 0; ti < 256; ti++) {\n"
			"\t\t\terr = 0.0;\n"
			"\t\t\tmax = 0.0;\n"
			"\t\t\tstd::vector<std::complex<double>> X(N/2+1);\n"
			"\t\t\tstd::vector<std::complex<double>> Y(N/2+1);\n"
			"\t\t\tstd::vector<double> y(N);\n"
			"\t\t\tstd::vector<double> x(N);\n"
			"\t\t\tfor (int n = 0; n < N / 2 + 1; n++) {\n"
			"\t\t\t\tX[n] = std::complex<double>(rand1(), rand1());\n"
			"\t\t\t}\n"
			"\t\t\tY = X;\n"
			"\t\t\ttm1.start();\n"
			"\t\t\ttm3.start();\n"
			"\t\t\tFFT_real_inv(X.data(), x.data(), N);\n"
			"\t\t\ttm1.stop();\n"
			"\t\t\ttm3.stop();\n"
			"\t\t\ttm2.start();\n"
			"\t\t\ttm4.start();\n"
			"\t\t\tfftw_real_inv(Y, y);\n"
			"\t\t\ttm2.stop();\n"
			"\t\t\ttm4.stop();\n"
			"\t\t\tfor (int i = 0; i < x.size(); i++) {\n"
				"\t\t\t\ty[i] -= x[i];\n"
			"\t\t\t}\n"
			"\t\t\tfor (int n = 0; n < N; n++) {\n"
			"\t\t\t\terr += std::abs(y[n]) * std::abs(y[n]);\n"
			"\t\t\t\t//printf(\"%i %e %e\\n\", n, x[n], y[n]);\n"
			"\t\t\t\tmax = std::max(max, std::abs(x[n]));\n"
			"\t\t\t}\n"
			"\t\t\terr = sqrt(err / N) / max;\n"
			"\t\t}\n"
			"\t\tprintf(\"%4i %4i %e %e %e %e %e %e %e\\n\", N, fft_nops_real[N], err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());\n"
			"");
	deindent();
	print("}\n");
	print("for( int N = 2; N <= %i; N+=%i) {\n", MAXFFT, DFFT);
	indent();
	print("if( !allow_real(N) ) {\n");
	indent();
	print("continue;\n");
	deindent();
	print("}\n");
	print_notab("\t\ttimer tm1, tm2;\n"
			"\t\tdouble err;\n"
			"\t\tdouble max;\n"
			"\t\tfor (int ti = 0; ti < 256; ti++) {\n"
			"\t\t\terr = 0.0;\n"
			"\t\t\tmax = 0.0;\n"
			"\t\t\tstd::vector<std::complex<double>> X(N/2+1);\n"
			"\t\t\tstd::vector<std::complex<double>> Y(N/2+1);\n"
			"\t\t\tstd::vector<double> y(N);\n"
			"\t\t\tstd::vector<double> x(N);\n"
			"\t\t\tfor (int n = 0; n < N; n++) {\n"
			"\t\t\t\tx[n] = rand1();\n"
			"\t\t\t}\n"
			"\t\t\ty = x;\n"
			"\t\t\ttm1.start();\n"
			"\t\t\ttm3.start();\n"
			"\t\t\tFFT_real(X.data(), x.data(), N);\n"
			"\t\t\ttm1.stop();\n"
			"\t\t\ttm3.stop();\n"
			"\t\t\ttm2.start();\n"
			"\t\t\ttm4.start();\n"
			"\t\t\tfftw_real(Y, y);\n"
			"\t\t\ttm2.stop();\n"
			"\t\t\ttm4.stop();\n"
			"\t\t\tfor (int i = 0; i < X.size(); i++) {\n"
			"\t\t\t\tY[i] -= X[i];\n"
			"\t\t\t}\n"
			"\t\t\tfor (int n = 0; n < N / 2 + 1; n++) {\n"
			"\t\t\t\terr += std::abs(Y[n]) * std::abs(Y[n]);\n"
			"\t\t\t\t//printf(\"%%i %%e %%e %%e %%e\\n\", n, X[n].real(), X[n].imag(), Y[n].real(), Y[n].imag());\n"
			"\t\t\t\tmax = std::max(max, std::abs(X[n]));\n"
			"\t\t\t}\n"
			"\t\t\terr = sqrt(err / N) / max;\n"
			"\t\t}\n"
			"\t\tprintf(\"%4i %4i %e %e %e %e %e %e %e\\n\", N, fft_nops_real[N], err, tm1.read(), tm2.read(), tm1.read() / tm2.read(), tm3.read(), tm4.read(), tm3.read() / tm4.read());\n"
			"");
	deindent();
	print("}\n");
	deindent();
	print("}\n");

	print("int main() {\n");
	indent();
	print("feenableexcept(FE_DIVBYZERO);\n");
	print("feenableexcept(FE_OVERFLOW);\n");
	print("feenableexcept(FE_INVALID);\n");
	print("test();\n");
	print("printf(\"\\n\");\n");
	print("test();\n");
	print("printf(\"EXITING\\n\");\n");
deindent();
	print("}\n\n");

	set_file("Makefile");
	print("CC=g++\n");
//	print("CFLAGS=-I. -g -O0 -D_GLIBCXX_DEBUG -fsanitize=address -static-libasan -march=native\n");
	print("CFLAGS=-I. -Ofast -march=native\n");
	print("DEPS = fft.hpp\n");
	print("OBJ = fft.o ");
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print_notab("fft.%i.o ", n);
	}
	for (int n = 2; n <= MAXFFT; n++) {
		print_notab("fft.real.%i.o ", n);
	}
	for (int n = 2; n <= MAXFFT; n++) {
		print_notab("fft.real.inv.%i.o ", n);
	}

	print("\n%.o: %.cpp $(DEPS)\n");
	print("\t$(CC) -c -o $@ $< $(CFLAGS)\n\n");

	print("ffttest: $(OBJ)\n");
	print("\t$(CC) -o $@ $^ $(CFLAGS) -lfftw3\n");
	set_file("");
	return 0;
}
