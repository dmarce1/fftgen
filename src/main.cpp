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
	print("\nvoid fft_base_%i(double* x) {\n", N);
	indent();
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	fft(N, 0, true);
	deindent();
	print("}\n\n");
	print("\nvoid fft_bitreverse_%i(double* x) {\n", N);
	indent();
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	fft_bitreverse(N);
	deindent();
	print("}\n\n");
	print("\nvoid fft_%i(double* x) {\n", N);
	indent();
	print("fft_bitreverse_%i(x);\n", N);
	print("fft_base_%i(x);\n", N);
	deindent();
	print("}\n\n");
	fft_nops[N] = fft_opcnt(N, 0);
}

void print_fft_real(int N) {
	std::string fname = "fft.real." + std::to_string(N) + ".cpp";
	set_file(fname);
	print("#include \"fft.hpp\"\n");
	print("\nvoid fft_real_base_%i(double* x) {\n", N);
	indent();
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	fft_real(N);
	deindent();
	print("}\n\n");

	print("\nvoid fft_real_%i(double* x, double* y) {\n", N);
	indent();
	for (int n = 0; n < NPAR; n++) {
		print("double tmp%i;\n", n);
	}
	fft_bitreverse_real(N);
	print("fft_real_base_%i(x);\n", N);
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
	print("\nvoid FFT_real_inv(std::complex<double>* xout, double* xin, int N);\n");
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
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print("void fft_base_%i(double*);\n", n);
	}
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print("void fft_real_base_%i(double*);\n", n);
	}
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print("void fft_%i(double*);\n", n);
	}
	for (int n = 2; n <= MAXFFT; n++) {
		print("void fft_real_%i(double*, double*);\n", n);
	}
	print("\n");

	set_file("fft.cpp");
	printf( "Hello World\n");
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

	print("typedef void (*func_real_type)(double*, double*);\n");
	print("typedef void (*func_type)(double*);\n");
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
	print("#define MAXFFT %i\n", MAXFFT);
	include("../gen_src/fft.cpp");

	set_file("Makefile");
	print("CC=g++\n");
//	print("CFLAGS=-I. -g -O0 -D_GLIBCXX_DEBUG -march=native\n");
	print("CFLAGS=-I. -Ofast -march=native\n");
	print("DEPS = fft.hpp\n");
	print("OBJ = fft.o ");
	for (int n = 2; n <= MAXFFT; n += DFFT) {
		print_notab("fft.%i.o ", n);
	}
	for (int n = 2; n <= MAXFFT; n++) {
		print_notab("fft.real.%i.o ", n);
	}

	print("\n%.o: %.cpp $(DEPS)\n");
	print("\t$(CC) -c -o $@ $< $(CFLAGS)\n\n");

	print("ffttest: $(OBJ)\n");
	print("\t$(CC) -o $@ $^ $(CFLAGS) -lfftw3\n");
	set_file("");
	return 0;
}
