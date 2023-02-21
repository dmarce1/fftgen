#include "fftgen.hpp"
#include <list>

void ct_fft(int N1, int N2, int o) {
	printf("// cooley-tukey - %i = %i x %i\n", N1 * N2, N1, N2);
	print("{\n");
	indent();
	int N = N1 * N2;
	for (int n1 = 0; n1 < N1; n1++) {
		fft(N2, n1 * N2 + o);
	}
	print("std::array<std::array<double, %i>, %i> z;\n", 2 * N1, N2);
	auto I = fft_bitreverse_indices(N1);
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			const auto W = twiddle(k2 * I[n1], N);
			if (k2 * I[n1] == 0) {
				print("z[%i][%i] = x[%i];\n", k2, 2 * n1 + 0, index(N2 * I[n1] + o, k2, 0, N2));
				print("z[%i][%i] = x[%i];\n", k2, 2 * n1 + 1, index(N2 * I[n1] + o, k2, 1, N2));
			} else if (k2 * I[n1] == N / 8 && N % 8 == 0) {
				print("z[%i][%i] = M_SQRT1_2 * (x[%i] + x[%i]);\n", k2, 2 * n1 + 0, index(N2 * I[n1] + o, k2, 0, N2), index(N2 * I[n1] + o, k2, 1, N2));
				print("z[%i][%i] = -M_SQRT1_2 * (x[%i] - x[%i]);\n", k2, 2 * n1 + 1, index(N2 * I[n1] + o, k2, 0, N2), index(N2 * I[n1] + o, k2, 1, N2));
			} else if (k2 * I[n1] == 3 * N / 8 && N % 8 == 0) {
				print("z[%i][%i] = -M_SQRT1_2 * (x[%i] - x[%i]);\n", k2, 2 * n1 + 0, index(N2 * I[n1] + o, k2, 0, N2), index(N2 * I[n1] + o, k2, 1, N2));
				print("z[%i][%i] = -M_SQRT1_2 * (x[%i] + x[%i]);\n", k2, 2 * n1 + 1, index(N2 * I[n1] + o, k2, 0, N2), index(N2 * I[n1] + o, k2, 1, N2));
			} else if (k2 * I[n1] == 5 * N / 8 && N % 8 == 0) {
				print("z[%i][%i] = -M_SQRT1_2 * (x[%i] + x[%i]);\n", k2, 2 * n1 + 0, index(N2 * I[n1] + o, k2, 0, N2), index(N2 * I[n1] + o, k2, 1, N2));
				print("z[%i][%i] = M_SQRT1_2 * (x[%i] - x[%i]);\n", k2, 2 * n1 + 1, index(N2 * I[n1] + o, k2, 0, N2), index(N2 * I[n1] + o, k2, 1, N2));
			} else if (k2 * I[n1] == 7 * N / 8 && N % 8 == 0) {
				print("z[%i][%i] = M_SQRT1_2 * (x[%i] - x[%i]);\n", k2, 2 * n1 + 0, index(N2 * I[n1] + o, k2, 0, N2), index(N2 * I[n1] + o, k2, 1, N2));
				print("z[%i][%i] = M_SQRT1_2 * (x[%i] + x[%i]);\n", k2, 2 * n1 + 1, index(N2 * I[n1] + o, k2, 0, N2), index(N2 * I[n1] + o, k2, 1, N2));
			} else if (k2 * I[n1] == N / 2 && N % 2 == 0) {
				print("z[%i][%i] = -x[%i];\n", k2, 2 * n1 + 0, index(N2 * I[n1] + o, k2, 0, N2));
				print("z[%i][%i] = -x[%i];\n", k2, 2 * n1 + 1, index(N2 * I[n1] + o, k2, 1, N2));
			} else if (k2 * I[n1] == N / 4 && N % 4 == 0) {
				print("z[%i][%i] = x[%i];\n", k2, 2 * n1 + 0, index(N2 * I[n1] + o, k2, 1, N2));
				print("z[%i][%i] = -x[%i];\n", k2, 2 * n1 + 1, index(N2 * I[n1] + o, k2, 0, N2));
			} else if (k2 * I[n1] == 3 * N / 4 && N % 4 == 0) {
				print("z[%i][%i] = -x[%i];\n", k2, 2 * n1 + 0, index(N2 * I[n1] + o, k2, 1, N2));
				print("z[%i][%i] = x[%i];\n", k2, 2 * n1 + 1, index(N2 * I[n1] + o, k2, 0, N2));
			} else {
				print("z[%i][%i] = std::fma(x[%i], (%.17e), x[%i] * (%.17e));\n", k2, 2 * n1 + 0, index(N2 * I[n1] + o, k2, 0, N2), W.real(), index(N2 * I[n1] + o, k2, 1, N2), -W.imag());
				print("z[%i][%i] = std::fma(x[%i], (%.17e), x[%i] * (%.17e));\n", k2, 2 * n1 + 1, index(N2 * I[n1] + o, k2, 0, N2), W.imag(), index(N2 * I[n1] + o, k2, 1, N2), W.real());
			}
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		print("{\n");
		indent();
		print("auto* x = z[%i].data();\n", k2);
		fft(N1, 0);
		deindent();
		print("}\n");
	}
	for (int k2 = 0; k2 < N2; k2++) {
		for (int k1 = 0; k1 < N1; k1++) {
			print("x[%i] = z[%i][%i];\n", index(o + N2 * k1, k2, 0, N), k2, 2 * k1 + 0);
			print("x[%i] = z[%i][%i];\n", index(o + N2 * k1, k2, 1, N), k2, 2 * k1 + 1);
		}
	}
	deindent();
	print("}\n");
}

int ct_fft_opcnt(int N1, int N2) {
	int cnt = 0;
	int N = N1 * N2;
	for (int n1 = 0; n1 < N1; n1++) {
		cnt += fft_opcnt(N2, 0);
	}
	auto I = fft_bitreverse_indices(N1);
	for (int k2 = 0; k2 < N2; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			const auto W = twiddle(k2 * I[n1], N);
			if (k2 * I[n1] == 0) {
			} else if (k2 * I[n1] == N / 2 && N % 2 == 0) {
				cnt += 2;
			} else if (k2 * I[n1] == N / 4 && N % 4 == 0) {
				cnt += 1;
			} else if (k2 * I[n1] == 3 * N / 4 && N % 4 == 0) {
				cnt += 1;
			} else {
				cnt += 4;
			}
		}
	}
	for (int k2 = 0; k2 < N2; k2++) {
		cnt += fft_opcnt(N1, 0);
	}
	return cnt;
}

std::vector<int> ct_fft_bitr(int N1, int N2, int o, std::vector<int> I) {
	printf("--?? %i %i\n", N1, N2);
	if (N1 == 8 && N2 == 14) {
		abort();
	}
	int r = N1;
	int N = N1 * N2;
	std::vector<int> L;
	for (int n = 0; n < r; n++) {
		std::vector<int> J;
		for (int k = 0; k < N / r; k++) {
			J.push_back(I[n + r * k]);
		}
		auto K = fft_bitr(N / r, o + n * N / r, J);
		L.insert(L.end(), K.begin(), K.end());
	}
	return L;
}

void print_z(int zi, int twi, int loci, int k, int r, int N, int o, bool ref) {
	const auto W = twiddle(twi * k, N);
	const int i = k + loci * N / r;
	/*if (twi * k == 0) {*/
	print("const auto%s zr%i = x[%i];\n", ref ? "&" : "", zi, index(o, i, 0, N));
	print("const auto%s zi%i = x[%i];\n", ref ? "&" : "", zi, index(o, i, 1, N));
	/*	} else if (twi * k == N / 8 && N % 8 == 0) {
	 print("const auto zr%i = M_SQRT1_2 * (x[%i] + x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	 print("const auto zi%i = -M_SQRT1_2 * (x[%i] - x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	 } else if (twi * k == N / 4 && N % 4 == 0) {
	 print("const auto%s zr%i = x[%i];\n", ref ? "&" : "", zi, index(o, i, 1, N));
	 print("const auto zi%i = -x[%i];\n", zi, index(o, i, 0, N));
	 } else if (twi * k == 3 * N / 8 && N % 8 == 0) {
	 print("const auto zr%i = -M_SQRT1_2 * (x[%i] - x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	 print("const auto zi%i = -M_SQRT1_2 * (x[%i] + x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	 } else if (twi * k == N / 2 && N % 2 == 0) {
	 print("const auto zr%i = -x[%i];\n", zi, index(o, i, 0, N));
	 print("const auto zi%i = -x[%i];\n", zi, index(o, i, 1, N));
	 } else if (twi * k == 5 * N / 8 && N % 8 == 0) {
	 print("const auto zr%i = -M_SQRT1_2 * (x[%i] + x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	 print("const auto zi%i = M_SQRT1_2 * (x[%i] - x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	 } else if (twi * k == 3 * N / 4 && N % 4 == 0) {
	 print("const auto zr%i = -x[%i];\n", zi, index(o, i, 1, N));
	 print("const auto%s zi%i = x[%i];\n", ref ? "&" : "", zi, index(o, i, 0, N));
	 } else if (twi * k == 7 * N / 8 && N % 8 == 0) {
	 print("const auto zr%i = M_SQRT1_2 * (x[%i] - x[%i];\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	 print("const auto zi%i = M_SQRT1_2 * (x[%i] + x[%i]);\n", zi, index(o, i, 0, N), index(o, i, 1, N));
	 } else {
	 print("const auto zr%i = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", zi, W.real(), index(o, i, 0, N), -W.imag(), index(o, i, 1, N));
	 print("const auto zi%i = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", zi, W.real(), index(o, i, 1, N), W.imag(), index(o, i, 0, N));
	 }*/
}

int print_z_opcnt(int zi, int twi, int loci, int k, int r, int N, int o) {
	int cnt = 0;
	const auto W = twiddle(twi * k, N);
	const int i = k + loci * N / r;
	if (twi * k == 0) {
	} else if (twi * k == N / 4 && N % 4 == 0) {
		cnt++;
	} else if (twi * k == N / 2 && N % 2 == 0) {
		cnt += 2;
	} else if (twi * k == 3 * N / 4 && N % 4 == 0) {
		cnt++;
	} else {
		cnt += 4;
	}
	return cnt;
}

std::vector<int> reverse_indices(std::vector<int> in) {
	int N = in.size();
	std::vector<int> out(N);
	for (int n = 0; n < N; n++) {
		out[in[n]] = n;
	}
	return out;
}

void fft_bitreverse(int N, std::vector<int> indices, int o) {
	if (indices.size() == 0) {
		indices = fft_bitreverse_indices(N);
//		indices = reverse_indices(indices);
	}
	std::vector<std::list<int>> strings;
	for (int l = 0; l < 2; l++) {
		std::vector<bool> touched(N, false);
		for (int k = 0; k < N; k++) {
			if (!touched[k]) {
				std::list<int> string;
				int current = k;
				int first = current;
				int next = indices[current];
				touched[current] = true;
				if (first != next) {
					string.push_back(-1);
					string.push_back(index(o, current, l, N));
					while (first != next) {
						string.push_back(index(o, next, l, N));
						current = next;
						touched[current] = true;
						next = indices[next];
					}
				}
				strings.push_back(string);
			}
		}
	}
	int n = 0;
	std::vector<std::vector<std::string>> cmds(NPAR);
	while (n < strings.size()) {
		int smallest = 9999999999, smalli;
		for (int i = 0; i < NPAR; i++) {
			if (cmds[i].size() < smallest) {
				smallest = cmds[i].size();
				smalli = i;
			}
		}
		std::vector<std::string> cmd;
		char* buf;
		while (strings[n].size()) {
			if (strings[n].front() == -1) {
				strings[n].pop_front();
				if (strings[n].size()) {
					asprintf(&buf, "tmp%i = x[%i];\n", smalli, strings[n].front());
					cmd.push_back(buf);
					free(buf);
				}
			} else {
				int to = strings[n].front();
				strings[n].pop_front();
				if (strings[n].size()) {
					asprintf(&buf, "x[%i] = x[%i];\n", to, strings[n].front());
					cmd.push_back(buf);
					free(buf);
				} else {
					asprintf(&buf, "x[%i] = tmp%i;\n", to, smalli);
					cmd.push_back(buf);
					free(buf);
				}
			}
		}
		cmds[smalli].insert(cmds[smalli].end(), cmd.begin(), cmd.end());
		n++;
	}
	n = 0;
	bool done = false;
	while (!done) {
		done = true;
		for (int j = 0; j < NPAR; j++) {
			if (cmds[j].size() > n) {
				print("%s", cmds[j][n].c_str());
				done = false;
			}
		}
		n++;
	}
}

std::vector<int> fft_radix_bitr(int r, int N, int o, std::vector<int> I) {
	std::vector<int> L;
	std::vector<int> J;
	std::vector<int> K;
	switch (r) {
	case 1:
		return I;
	default:
		for (int n = 0; n < r; n++) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + r * k]);
			}
			auto K = fft_bitr(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
	};
	return L;
}

std::vector<int> fft_bitr(int N, int o, std::vector<int> indices, bool first) {
	if (N == 1) {
		return indices;
	}
	auto fftt = best_radix(N, o, first);
	if (fftt.type == RADIX) {
		return fft_radix_bitr(fftt.N1, N, o, indices);
	} else if (fftt.type == GOOD) {
		return gt2_fft_bitr(fftt.N1, fftt.N2, o, indices);
	} else if (fftt.type == COOLEY) {
		printf("----%i %i\n", fftt.N1, fftt.N2);
		return ct_fft_bitr(fftt.N1, fftt.N2, o, indices);
	} else {
		std::vector<int> I3(N);
		std::vector<int> I1(N);
		const auto gq = raders_gq(N);
		auto I2 = fft_bitreverse_indices(N - 1);
		for (int i = 0; i < I2.size(); i++) {
			I3[i] = indices[gq[I2[i]]];
		}
		I3[N - 1] = indices[0];
		return I3;
	}
}

std::vector<int> fft_bitreverse_indices(int N) {
	std::vector<int> I;
	std::vector<int> indices;
	I.resize(0);
	for (int n = 0; n < N; n++) {
		I.push_back(n);
	}
	indices = fft_bitr(N, 0, I, true);
	return indices;
}
void fft_radix_dit(int r, int N, int o);
void fft_radix_dif(int r, int N, int o);

void fft_radix(int r, int N, int o) {
	return fft_radix_dit(r, N, o);
}

void fft_radix_dit(int r, int N, int o) {
	print("{\n");
	print("// radix - %i\n", r);
	indent();
	if (N > r) {
		print("for (int n = 0; n < %i; n++) {\n", r);
		indent();
		print("fft_base_%i(x + (%i * n) / %i);\n", N / r, 2 * N, r);
		deindent();
		print("}\n");
	}
	if (N >= LOOP_N) {
		if (N / r > 1) {
			print("static twiddle_set twiddles(%i);\n", N);
			print("for( int n1 = 1; n1 < %i; n1++ ) {\n", r);
			indent();
			print("const int n1N2 = n1 * %i;\n", N / r);
			print("for( int k2 = 1; k2 < %i; k2++ ) {\n", N / r);
			indent();
			print("const int i = n1N2 + k2;\n");
			print("const auto tw = twiddles[n1 * k2];\n");
			print("const int ir = 2 * i;\n");
			print("const int ii = ir + 1;\n");
			print("tmp0 = x[ir];\n");
			print("x[ir] = x[ir] * tw.real() - x[ii] * tw.imag();\n");
			print("x[ii] = std::fma(tmp0, tw.imag(), x[ii] * tw.real());\n");
			deindent();
			print("}\n");
			deindent();
			print("}\n");
		}
	} else {
		for (int i = 0; i < N; i++) {
			const int k = i % (N / r);
			const int n = i / (N / r);
			int ir, ii;
			ir = index(o, k + n * N / r, 0, N / r);
			ii = index(o, k + n * N / r, 1, N / r);
			int tmpi = i % NPAR;
			if (k * n == 0) {
			} else if (k * n == N / 4 && N % 4 == 0) {
				print("std::swap(x[%i], x[%i]);\n", ir, ii);
				print("x[%i] = -x[%i];\n", ii, ii);
			} else if (k * n == 3 * N / 4 && N % 4 == 0) {
				print("std::swap(x[%i], x[%i]);\n", ir, ii);
				print("x[%i] = -x[%i];\n", ir, ir);
			} else if (k * n == N / 2 && N % 2 == 0) {
				print("x[%i] = -x[%i];\n", ir, ir);
				print("x[%i] = -x[%i];\n", ii, ii);
			} else {
				const auto W = twiddle(k * n, N);
				print("tmp%i = x[%i];\n", tmpi, ir);
				print("x[%i] = std::fma(x[%i], (%.17e), x[%i] * (%.17e));\n", ir, ir, W.real(), ii, -W.imag());
				print("x[%i] = std::fma(x[%i], (%.17e), tmp%i * (%.17e));\n", ii, ii, W.real(), tmpi, W.imag());
			}
		}
	}

	if (N / r != 1) {
		print("for (int k = 0; k < %i; k++) {\n", N / r);
		indent();
		print("auto* y = x + 2 * k;\n");
	} else {
		print("{\n");
		indent();
		print("auto* y = x;\n");
	}
	std::vector<std::string> out;
	std::vector<std::string> in;
	for (int i = 0; i < r; i++) {
		out.push_back(std::string("y[" + std::to_string(2 * (i * N / r)) + "]"));
		out.push_back(std::string("y[" + std::to_string(2 * (i * N / r) + 1) + "]"));
	}
	in = out;
	print_complex_short_fft(r, in, out);
	deindent();
	print("}\n");
	deindent();
	print("}\n");
}

void fft_radix_dit_strided(int r, int N, int o, int s) {
	print("{\n");
	print("// radix - %i\n", r);
	indent();
	if (N > r) {
		int R = best_radix(N / r, 0).N1;
		for (int n = 0; n < r; n++) {
			fft_radix_dit_strided(R, N / r, n * N / r, s);
		}
	}
	if (N >= LOOP_N) {
		if (N / r > 1) {
			print("static twiddle_set twiddles(%i);\n", N);
			print("for( int n1 = 1; n1 < %i; n1++ ) {\n", r);
			indent();
			print("const int n1N2 = n1 * %i;\n", N / r);
			print("for( int k2 = 1; k2 < %i; k2++ ) {\n", N / r);
			indent();
			print("const int i = n1N2 + k2;\n");
			print("const auto tw = twiddles[n1 * k2];\n");
			print("const int ir = %i * i;\n", 2 * s + 2 * o);
			print("const int ii = ir + 1;\n");
			print("tmp0 = x[ir];\n");
			print("x[ir] = x[ir] * tw.real() - x[ii] * tw.imag();\n");
			print("x[ii] = std::fma(tmp0, tw.imag(), x[ii] * tw.real());\n");
			deindent();
			print("}\n");
			deindent();
			print("}\n");
		}
	} else {
		for (int i = 0; i < N; i++) {
			const int k = i % (N / r);
			const int n = i / (N / r);
			int ir, ii;
			ir = index(o, s * (k + n * N / r), 0, N / r);
			ii = index(o, s * (k + n * N / r), 1, N / r);
			int tmpi = i % NPAR;
			if (k * n == 0) {
			} else if (k * n == N / 8 && N % 8 == 0) {
				print("tmp%i = x[%i];\n", tmpi, ir);
				print("x[%i] = M_SQRT1_2 * (x[%i] + x[%i]);\n", ir, ir, ii);
				print("x[%i] = -M_SQRT1_2 * (tmp%i - x[%i]);\n", ii, tmpi, ii);
			} else if (k * n == 3 * N / 8 && N % 8 == 0) {
				print("tmp%i = x[%i];\n", tmpi, ir);
				print("x[%i] = -M_SQRT1_2 * (x[%i] - x[%i]);\n", ir, ir, ii);
				print("x[%i] = -M_SQRT1_2 * (tmp%i + x[%i]);\n", ii, tmpi, ii);
			} else if (k * n == 7 * N / 8 && N % 8 == 0) {
				print("tmp%i = x[%i];\n", tmpi, ir);
				print("x[%i] = M_SQRT1_2 * (x[%i] - x[%i]);\n", ir, ir, ii);
				print("x[%i] = M_SQRT1_2 * (tmp%i + x[%i]);\n", ii, tmpi, ii);
			} else if (k * n == 5 * N / 8 && N % 8 == 0) {
				print("tmp%i = x[%i];\n", tmpi, ir);
				print("x[%i] = -M_SQRT1_2 * (x[%i] + x[%i]);\n", ir, ir, ii);
				print("x[%i] = M_SQRT1_2 * (tmp%i - x[%i]);\n", ii, tmpi, ii);
			} else if (k * n == N / 4 && N % 4 == 0) {
				print("std::swap(x[%i], x[%i]);\n", ir, ii);
				print("x[%i] = -x[%i];\n", ii, ii);
			} else if (k * n == 3 * N / 4 && N % 4 == 0) {
				print("std::swap(x[%i], x[%i]);\n", ir, ii);
				print("x[%i] = -x[%i];\n", ir, ir);
			} else if (k * n == N / 2 && N % 2 == 0) {
				print("x[%i] = -x[%i];\n", ir, ir);
				print("x[%i] = -x[%i];\n", ii, ii);
			} else {
				const auto W = twiddle(k * n, N);
				print("tmp%i = x[%i];\n", tmpi, ir);
				print("x[%i] = std::fma(x[%i], (%.17e), x[%i] * (%.17e));\n", ir, ir, W.real(), ii, -W.imag());
				print("x[%i] = std::fma(x[%i], (%.17e), tmp%i * (%.17e));\n", ii, ii, W.real(), tmpi, W.imag());
			}
		}
	}

	if (N / r != 1) {
		print("for (int k = 0; k < %i; k++) {\n", N / r);
		indent();
		print("auto* y = x + %i * k + %i;\n", s * 2, 2 * o);
	} else {
		print("{\n");
		indent();
		print("auto* y = x + %i;\n", 2 * o);
	}
	std::vector<std::string> out;
	std::vector<std::string> in;
	for (int i = 0; i < r; i++) {
		out.push_back(std::string("y[" + std::to_string((s * 2 * i * N / r)) + "]"));
		out.push_back(std::string("y[" + std::to_string((s * 2 * i * N / r) + 1) + "]"));
	}
	in = out;
	print_complex_short_fft(r, in, out);
	deindent();
	print("}\n");
	deindent();
	print("}\n");
}



fft_type best_radix(int N, int o, bool first) {

	fft_type fftt;

	auto facs = fft_factorization(N);
	int radix = facs[0].first;
	fftt.nops = 0;
	if( radix > RADER_LEN) {
		fftt.type = RADERS;
		fftt.N1 = radix;
	} else {
		fftt.type = RADIX;
		fftt.N1 = radix;
	}
	return fftt;
}

int fft_opcnt(int N, bool first) {
	auto i = best_radix(N, 0, first);
	return i.nops;
}

void fft(int N, int o, bool first) {
	if (N == 1) {
		return;
	}
	auto fftt = best_radix(N, o, first);
	if (fftt.type == RADIX) {
		fft_radix(fftt.N1, N, o);
	}
	if (fftt.type == RADERS) {
		raders_fft(N, o);
	}
	if (fftt.type == GOOD) {
		gt2_fft(fftt.N1, fftt.N2, o);
	}
	if (fftt.type == COOLEY) {
		ct_fft(fftt.N1, fftt.N2, o);
	}
}
int fft_radix_opcnt(int r, int N) {

	if (N < r || r == 1) {
		return 0;
	}

	int cnt = 0;
	if (N > r) {
		if (r == 4) {
			cnt += fft_opcnt(N / 2);
			cnt += 2 * fft_opcnt(N / 4);
		} else {
			for (int i = 0; i < r; i++) {
				cnt += fft_opcnt(N / r);
			}
		}
	}
	for (int k = 0; k < N / r; k++) {
		if (r == 4) {
			cnt += print_z_opcnt(1, 1, 2, k, r, N, 0);
			cnt += print_z_opcnt(3, 3, 3, k, r, N, 0);
		} else {
			for (int i = 0; i < r; i++) {
				cnt += print_z_opcnt(i, i, i, k, r, N, 0);
			}
		}
	}
	switch (r) {
	case 2:
		cnt += N / r * 4;
		return cnt;
	case 3:
		cnt += N / r * 14;
		return cnt;
	case 4:
		cnt += N / r * 12;
		return cnt;
	case 5:
		cnt += N / r * 34;
		return cnt;
	case 6:
		cnt += N / r * 40;
		return cnt;
	case 7:
		cnt += N / r * 40;
		return cnt;
	case 8:
		cnt += N / r * 40;
		return cnt;
	default:
		for (int k = 0; k < N / r; k++) {
			for (int i = 1; i <= (r - 1) / 2; i++) {
				cnt += 4;
				for (int j = 1; j <= (r - 1) / 2; j++) {
					cnt += 4;
				}
			}
			for (int i = 1; i < r; i++) {
				cnt += 2;
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				cnt += 4;
			}
		}
		return cnt;
	};
}
