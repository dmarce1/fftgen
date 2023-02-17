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
		for (int n = 0; n < r; n++) {
			print("fft_base_%i(x + %i);\n", N / r, 2 * (o + n * N / r));
		}
	}
	if (N >= LOOP_N) {
		if (N / r > 1) {
			print("static twiddle_set<%i> twiddles;\n", N);
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
		print("auto* y = x + 2 * k;\n");
	} else {
		print("{\n");
		indent();
		print("auto* y = x;\n");
	}
	for (int i = 0; i < r; i++) {
		print("const auto zr%i = y[%i];\n", i, index(o, i * N / r, 0, N / r));
		print("const auto zi%i = y[%i];\n", i, index(o, i * N / r, 1, N / r));
	}

	switch (r) {
	case 1:
		return;
	case 2:
		print("y[%i] = zr0 + zr1;\n", index(o, 0, 0, N));
		print("y[%i] = zi0 + zi1;\n", index(o, 0, 1, N));
		print("y[%i] = zr0 - zr1;\n", index(o, 0 + N / 2, 0, N));
		print("y[%i] = zi0 - zi1;\n", index(o, 0 + N / 2, 1, N));
		break;
	case 3:
		print("const auto tr1 = zr1 + zr2;\n");
		print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
		print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
		print("const auto ti1 = zi1 + zi2;\n");
		print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
		print("const auto ti3 = (%24.17e) * (zi1 - zi2);\n", sqrt(3) * 0.5);
		print("y[%i] = zr0 + tr1;\n", index(o, 0, 0, N));
		print("y[%i] = zi0 + ti1;\n", index(o, 0, 1, N));
		print("y[%i] = tr2 + ti3;\n", index(o, 0 + N / 3, 0, N));
		print("y[%i] = ti2 - tr3;\n", index(o, 0 + N / 3, 1, N));
		print("y[%i] = tr2 - ti3;\n", index(o, 0 + 2 * N / 3, 0, N));
		print("y[%i] = ti2 + tr3;\n", index(o, 0 + 2 * N / 3, 1, N));
		break;
	case 4:
		print("const auto tr1 = zr0 + zr2;\n");
		print("const auto ti1 = zi0 + zi2;\n");
		print("const auto tr3 = zr0 - zr2;\n");
		print("const auto ti3 = zi0 - zi2;\n");
		print("const auto tr2 = zr1 + zr3;\n");
		print("const auto ti2 = zi1 + zi3;\n");
		print("const auto tr4 = zr1 - zr3;\n");
		print("const auto ti4 = zi1 - zi3;\n");
		print("y[%i] = tr1 + tr2;\n", index(o, 0, 0, N));
		print("y[%i] = ti1 + ti2;\n", index(o, 0, 1, N));
		print("y[%i] = tr3 + ti4;\n", index(o, 0 + N / 4, 0, N));
		print("y[%i] = ti3 - tr4;\n", index(o, 0 + N / 4, 1, N));
		print("y[%i] = tr1 - tr2;\n", index(o, 0 + N / 2, 0, N));
		print("y[%i] = ti1 - ti2;\n", index(o, 0 + N / 2, 1, N));
		print("y[%i] = tr3 - ti4;\n", index(o, 0 + 3 * N / 4, 0, N));
		print("y[%i] = ti3 + tr4;\n", index(o, 0 + 3 * N / 4, 1, N));
		break;
	case 5:
		print("const auto tr1 = zr1 + zr4;\n");
		print("const auto tr2 = zr2 + zr3;\n");
		print("const auto tr3 = zr1 - zr4;\n");
		print("const auto tr4 = zr2 - zr3;\n");
		print("const auto tr5 = tr1 + tr2;\n");
		print("const auto tr6 = (%24.17e) * (tr1 - tr2);\n", sqrt(5) * 0.25);
		print("const auto tr7 = std::fma(tr5, -0.25, zr0);\n");
		print("const auto tr8 = tr7 + tr6;\n");
		print("const auto tr9 = tr7 - tr6;\n");
		print("const auto tr10 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
		print("const auto tr11 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
		print("const auto ti1 = zi1 + zi4;\n");
		print("const auto ti2 = zi2 + zi3;\n");
		print("const auto ti3 = zi1 - zi4;\n");
		print("const auto ti4 = zi2 - zi3;\n");
		print("const auto ti5 = ti1 + ti2;\n");
		print("const auto ti6 = (%24.17e) * (ti1 - ti2);\n", sqrt(5) * 0.25);
		print("const auto ti7 = std::fma(ti5, -0.25, zi0);\n");
		print("const auto ti8 = ti7 + ti6;\n");
		print("const auto ti9 = ti7 - ti6;\n");
		print("const auto ti10 = std::fma((%24.17e), ti3, (%24.17e) * ti4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
		print("const auto ti11 = std::fma((%24.17e), ti3, (%24.17e) * ti4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
		print("y[%i] = zr0 + tr5;\n", index(o, 0 + 0 * N / 5, 0, N));
		print("y[%i] = zi0 + ti5;\n", index(o, 0 + 0 * N / 5, 1, N));
		print("y[%i] = tr8 + ti10;\n", index(o, 0 + 1 * N / 5, 0, N));
		print("y[%i] = ti8 - tr10;\n", index(o, 0 + 1 * N / 5, 1, N));
		print("y[%i] = tr9 + ti11;\n", index(o, 0 + 2 * N / 5, 0, N));
		print("y[%i] = ti9 - tr11;\n", index(o, 0 + 2 * N / 5, 1, N));
		print("y[%i] = tr9 - ti11;\n", index(o, 0 + 3 * N / 5, 0, N));
		print("y[%i] = ti9 + tr11;\n", index(o, 0 + 3 * N / 5, 1, N));
		print("y[%i] = tr8 - ti10;\n", index(o, 0 + 4 * N / 5, 0, N));
		print("y[%i] = ti8 + tr10;\n", index(o, 0 + 4 * N / 5, 1, N));
		break;
	case 6:
		print("const auto tr1 = zr2 + zr4;\n");
		print("const auto ti1 = zi2 + zi4;\n");
		print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
		print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
		print("const auto tr3 = (%24.17e) * (zr2 - zr4);\n", sin(M_PI / 3.0));
		print("const auto ti3 = (%24.17e) * (zi2 - zi4);\n", sin(M_PI / 3.0));
		print("const auto tr4 = zr5 + zr1;\n");
		print("const auto ti4 = zi5 + zi1;\n");
		print("const auto tr5 = std::fma(tr4, -0.5, zr3);\n");
		print("const auto ti5 = std::fma(ti4, -0.5, zi3);\n");
		print("const auto tr6 = (%24.17e) * (zr5 - zr1);\n", sin(M_PI / 3.0));
		print("const auto ti6 = (%24.17e) * (zi5 - zi1);\n", sin(M_PI / 3.0));
		print("const auto tr7 = zr0 + tr1;\n");
		print("const auto ti7 = zi0 + ti1;\n");
		print("const auto tr8 = tr2 + ti3;\n");
		print("const auto ti8 = ti2 - tr3;\n");
		print("const auto tr9 = tr2 - ti3;\n");
		print("const auto ti9 = ti2 + tr3;\n");
		print("const auto tr10 = zr3 + tr4;\n");
		print("const auto ti10 = zi3 + ti4;\n");
		print("const auto tr11 = tr5 + ti6;\n");
		print("const auto ti11 = ti5 - tr6;\n");
		print("const auto tr12 = tr5 - ti6;\n");
		print("const auto ti12 = ti5 + tr6;\n");
		print("y[%i] = tr7 + tr10;\n", index(o, 0 + 0 * N / 6, 0, N));
		print("y[%i] = tr8 - tr11;\n", index(o, 0 + 1 * N / 6, 0, N));
		print("y[%i] = tr9 + tr12;\n", index(o, 0 + 2 * N / 6, 0, N));
		print("y[%i] = tr7 - tr10;\n", index(o, 0 + 3 * N / 6, 0, N));
		print("y[%i] = tr8 + tr11;\n", index(o, 0 + 4 * N / 6, 0, N));
		print("y[%i] = tr9 - tr12;\n", index(o, 0 + 5 * N / 6, 0, N));
		print("y[%i] = ti7 + ti10;\n", index(o, 0 + 0 * N / 6, 1, N));
		print("y[%i] = ti8 - ti11;\n", index(o, 0 + 1 * N / 6, 1, N));
		print("y[%i] = ti9 + ti12;\n", index(o, 0 + 2 * N / 6, 1, N));
		print("y[%i] = ti7 - ti10;\n", index(o, 0 + 3 * N / 6, 1, N));
		print("y[%i] = ti8 + ti11;\n", index(o, 0 + 4 * N / 6, 1, N));
		print("y[%i] = ti9 - ti12;\n", index(o, 0 + 5 * N / 6, 1, N));
		break;
	case 7: {
		constexpr double theta = 2.0 * M_PI / 7.0;
		const double c1 = cos(theta);
		const double c2 = cos(2.0 * theta);
		const double c3 = cos(3.0 * theta);
		const double c4 = sin(theta);
		const double c5 = sin(2.0 * theta);
		const double c6 = sin(3.0 * theta);
		print("const auto tr1 = zr1 + zr6;\n");
		print("const auto ti1 = zi1 + zi6;\n");
		print("const auto tr2 = zr2 + zr5;\n");
		print("const auto ti2 = zi2 + zi5;\n");
		print("const auto tr3 = zr3 + zr4;\n");
		print("const auto ti3 = zi3 + zi4;\n");
		print("const auto tr4 = zr1 - zr6;\n");
		print("const auto ti4 = zi1 - zi6;\n");
		print("const auto tr5 = zr2 - zr5;\n");
		print("const auto ti5 = zi2 - zi5;\n");
		print("const auto tr6 = zr3 - zr4;\n");
		print("const auto ti6 = zi3 - zi4;\n");
		print("const auto tr7 = zr0 - 0.5 * tr3;\n");
		print("const auto ti7 = zi0 - 0.5 * ti3;\n");
		print("const auto tr8 = tr1 - tr3;\n");
		print("const auto ti8 = ti1 - ti3;\n");
		print("const auto tr9 = tr2 - tr3;\n");
		print("const auto ti9 = ti2 - ti3;\n");
		print("const auto yr0 = zr0 + tr1 + tr2 + tr3;\n");
		print("const auto yi0 = zi0 + ti1 + ti2 + ti3;\n");
		print("const auto yr1 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c1, c2);
		print("const auto yi1 = std::fma(ti8, %.17e, std::fma(ti9, %.17e, ti7));\n", c1, c2);
		print("const auto yr2 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c2, c3);
		print("const auto yi2 = std::fma(ti8, %.17e, std::fma(ti9, %.17e, ti7));\n", c2, c3);
		print("const auto yr3 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c3, c1);
		print("const auto yi3 = std::fma(ti8, %.17e, std::fma(ti9, %.17e, ti7));\n", c3, c1);
		print("const auto yr4 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c6, -c4, c5);
		print("const auto yi4 = std::fma(ti4, %.17e, std::fma(ti5, %.17e, ti6 * %.17e));\n", c6, -c4, c5);
		print("const auto yr5 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c5, -c6, -c4);
		print("const auto yi5 = std::fma(ti4, %.17e, std::fma(ti5, %.17e, ti6 * %.17e));\n", c5, -c6, -c4);
		print("const auto yr6 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c4, c5, c6);
		print("const auto yi6 = std::fma(ti4, %.17e, std::fma(ti5, %.17e, ti6 * %.17e));\n", c4, c5, c6);
		print("y[%i] = yr0;\n", index(o, 0 + 0 * N / 7, 0, N));
		print("y[%i] = yi0;\n", index(o, 0 + 0 * N / 7, 1, N));
		print("y[%i] = yr1 + yi6;\n", index(o, 0 + 1 * N / 7, 0, N));
		print("y[%i] = yi1 - yr6;\n", index(o, 0 + 1 * N / 7, 1, N));
		print("y[%i] = yr2 + yi5;\n", index(o, 0 + 2 * N / 7, 0, N));
		print("y[%i] = yi2 - yr5;\n", index(o, 0 + 2 * N / 7, 1, N));
		print("y[%i] = yr3 + yi4;\n", index(o, 0 + 3 * N / 7, 0, N));
		print("y[%i] = yi3 - yr4;\n", index(o, 0 + 3 * N / 7, 1, N));
		print("y[%i] = yr3 - yi4;\n", index(o, 0 + 4 * N / 7, 0, N));
		print("y[%i] = yi3 + yr4;\n", index(o, 0 + 4 * N / 7, 1, N));
		print("y[%i] = yr2 - yi5;\n", index(o, 0 + 5 * N / 7, 0, N));
		print("y[%i] = yi2 + yr5;\n", index(o, 0 + 5 * N / 7, 1, N));
		print("y[%i] = yr1 - yi6;\n", index(o, 0 + 6 * N / 7, 0, N));
		print("y[%i] = yi1 + yr6;\n", index(o, 0 + 6 * N / 7, 1, N));
	}
		break;
	case 8:
		printf("N = %i - real Radix 8\n", N);
		print("const auto tr1 = zr0 + zr4;\n");
		print("const auto ti1 = zi0 + zi4;\n");
		print("const auto tr2 = zr0 - zr4;\n");
		print("const auto ti2 = zi0 - zi4;\n");
		print("const auto tr3 = zr1 + zr5;\n");
		print("const auto ti3 = zi1 + zi5;\n");
		print("const auto tr4 = zr1 - zr5;\n");
		print("const auto ti4 = zi1 - zi5;\n");
		print("const auto tr5 = zr2 + zr6;\n");
		print("const auto ti5 = zi2 + zi6;\n");
		print("const auto tr6 = zr2 - zr6;\n");
		print("const auto ti6 = zi2 - zi6;\n");
		print("const auto tr7 = zr3 + zr7;\n");
		print("const auto ti7 = zi3 + zi7;\n");
		print("const auto tr8 = zr3 - zr7;\n");
		print("const auto ti8 = zi3 - zi7;\n");
		print("const auto tr9 = tr1 + tr5;\n");
		print("const auto ti9 = ti1 + ti5;\n");
		print("const auto tr10 = tr3 + tr7;\n");
		print("const auto ti10 = ti3 + ti7;\n");
		print("const auto tr11 = (%.17e) * (tr4 - tr8);\n", 1.0 / sqrt(2));
		print("const auto ti11 = (%.17e) * (ti4 - ti8);\n", 1.0 / sqrt(2));
		print("const auto tr12 = (%.17e) * (tr4 + tr8);\n", 1.0 / sqrt(2));
		print("const auto ti12 = (%.17e) * (ti4 + ti8);\n", 1.0 / sqrt(2));
		print("const auto yr0 = tr9 + tr10;\n");
		print("const auto yi0 = ti9 + ti10;\n");
		print("const auto yr1 = tr2 + tr11;\n");
		print("const auto yi1 = ti2 + ti11;\n");
		print("const auto yr2 = tr1 - tr5;\n");
		print("const auto yi2 = ti1 - ti5;\n");
		print("const auto yr3 = tr2 - tr11;\n");
		print("const auto yi3 = ti2 - ti11;\n");
		print("const auto yr4 = tr9 - tr10;\n");
		print("const auto yi4 = ti9 - ti10;\n");
		print("const auto yr5 = tr12 - tr6;\n");
		print("const auto yi5 = ti12 - ti6;\n");
		print("const auto yr6 = tr3 - tr7;\n");
		print("const auto yi6 = ti3 - ti7;\n");
		print("const auto yr7 = tr12 + tr6;\n");
		print("const auto yi7 = ti12 + ti6;\n");
		print("y[%i] = yr0;\n", index(o, 0 + 0 * N / 8, 0, N));
		print("y[%i] = yi0;\n", index(o, 0 + 0 * N / 8, 1, N));
		print("y[%i] = yr1 + yi7;\n", index(o, 0 + 1 * N / 8, 0, N));
		print("y[%i] = yi1 - yr7;\n", index(o, 0 + 1 * N / 8, 1, N));
		print("y[%i] = yr2 + yi6;\n", index(o, 0 + 2 * N / 8, 0, N));
		print("y[%i] = yi2 - yr6;\n", index(o, 0 + 2 * N / 8, 1, N));
		print("y[%i] = yr3 + yi5;\n", index(o, 0 + 3 * N / 8, 0, N));
		print("y[%i] = yi3 - yr5;\n", index(o, 0 + 3 * N / 8, 1, N));
		print("y[%i] = yr4;\n", index(o, 0 + 4 * N / 8, 0, N));
		print("y[%i] = yi4;\n", index(o, 0 + 4 * N / 8, 1, N));
		print("y[%i] = yr3 - yi5;\n", index(o, 0 + 5 * N / 8, 0, N));
		print("y[%i] = yi3 + yr5;\n", index(o, 0 + 5 * N / 8, 1, N));
		print("y[%i] = yr2 - yi6;\n", index(o, 0 + 6 * N / 8, 0, N));
		print("y[%i] = yi2 + yr6;\n", index(o, 0 + 6 * N / 8, 1, N));
		print("y[%i] = yr1 - yi7;\n", index(o, 0 + 7 * N / 8, 0, N));
		print("y[%i] = yi1 + yr7;\n", index(o, 0 + 7 * N / 8, 1, N));
		break;
		/*case 9:
		printf("Using radix 9 for n = %i\n", N);
		{
			const auto theta = 2.0 * M_PI / 9.0;
			const auto c1 = sin(3.0 * theta);
			const auto c2 = cos(theta);
			const auto c3 = sin(theta);
			const auto c4 = cos(2.0 * theta);
			const auto c5 = sin(2.0 * theta);
			const auto c6 = c1 * c2;
			const auto c7 = c1 * c3;
			const auto c8 = c1 * c4;
			const auto c9 = c1 * c5;
			print("const auto tr1 = zr3 + zr6;\n");
			print("y[%i] = yr0;\n", index(o, 0 + 0 * N / 9, 0, N));
			print("y[%i] = yi0;\n", index(o, 0 + 0 * N / 9, 1, N));
			print("y[%i] = yr1 + yi8;\n", index(o, 0 + 1 * N / 9, 0, N));
			print("y[%i] = yi1 - yr8;\n", index(o, 0 + 1 * N / 9, 1, N));
			print("y[%i] = yr2 + yi7;\n", index(o, 0 + 2 * N / 9, 0, N));
			print("y[%i] = yi2 - yr7;\n", index(o, 0 + 2 * N / 9, 1, N));
			print("y[%i] = yr3 + yi6;\n", index(o, 0 + 3 * N / 9, 0, N));
			print("y[%i] = yi3 - yr6;\n", index(o, 0 + 3 * N / 9, 1, N));
			print("y[%i] = yr4 + yi5;\n", index(o, 0 + 4 * N / 9, 0, N));
			print("y[%i] = yi4 - yr5;\n", index(o, 0 + 4 * N / 9, 1, N));
			print("y[%i] = yr4 - yi5;\n", index(o, 0 + 5 * N / 9, 0, N));
			print("y[%i] = yi4 + yr5;\n", index(o, 0 + 5 * N / 9, 1, N));
			print("y[%i] = yr3 - yi6;\n", index(o, 0 + 6 * N / 9, 0, N));
			print("y[%i] = yi3 + yr6;\n", index(o, 0 + 6 * N / 9, 1, N));
			print("y[%i] = yr2 - yi7;\n", index(o, 0 + 7 * N / 9, 0, N));
			print("y[%i] = yi2 + yr7;\n", index(o, 0 + 7 * N / 9, 1, N));
			print("y[%i] = yr1 - yi8;\n", index(o, 0 + 8 * N / 9, 0, N));
			print("y[%i] = yi1 + yr8;\n", index(o, 0 + 8 * N / 9, 1, N));
		}
		break;*/
	default:
		for (int j = 1; j <= (r - 1) / 2; j++) {
			print("const auto txp%i = zr%i + zr%i;\n", j, j, r - j);
		}
		for (int j = 1; j <= (r - 1) / 2; j++) {
			print("const auto txm%i = zr%i - zr%i;\n", j, j, r - j);
		}
		for (int j = 1; j <= (r - 1) / 2; j++) {
			print("const auto typ%i = zi%i + zi%i;\n", j, j, r - j);
		}
		for (int j = 1; j <= (r - 1) / 2; j++) {
			print("const auto tym%i = zi%i - zi%i;\n", j, j, r - j);
		}
		for (int i = 1; i <= (r - 1) / 2; i++) {
			print("auto ap%i = zr0;\n", i);
			print("auto bp%i = zi0;\n", i);
			for (int j = 1; j <= (r - 1) / 2; j++) {
				print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / r), i);
				print("bp%i = std::fma(typ%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / r), i);
				if (j == 1) {
					print("double am%i = tym%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
					print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
				} else {
					print("am%i = std::fma(tym%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
					print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
				}
			}
		}
		print("y[%i] = ", index(o, 0, 0, N));
		for (int i = 0; i < r; i++) {
			print_notab(" + zr%i", i);
		}
		print_notab(";\n");
		print("y[%i] = ", index(o, 0, 1, N));
		for (int i = 0; i < r; i++) {
			print_notab(" + zi%i", i);
		}
		print_notab(";\n");
		for (int i = 1; i <= (r - 1) / 2; i++) {
			print("y[%i] = ap%i + am%i;\n", index(o, 0 + i * N / r, 0, N), i, i);
			print("y[%i] = bp%i - bm%i;\n", index(o, 0 + i * N / r, 1, N), i, i);
			print("y[%i] = ap%i - am%i;\n", index(o, 0 + (r - i) * N / r, 0, N), i, i);
			print("y[%i] = bp%i + bm%i;\n", index(o, 0 + (r - i) * N / r, 1, N), i, i);
		}
	}
	deindent();
	print("}\n");
	deindent();
	print("}\n");
}

fft_type best_radix(int N, int o, bool first) {
	fft_type fftt;
	fftt.type = RADIX;
	fftt.nops = 0;
	const int rader_radix = is_prime(N) ? RADER_LEN : 1000000000;
	if (N % 6 == 0) {
		fftt.N1 = 6;
	} else {
		int twopow = 0;
		int n = N;
		while (n % 2 == 0) {
			twopow++;
			n /= 2;
		}
		if (twopow >= 3 && twopow % 3 == 0) {
			fftt.N1 = 8;
		} else if (twopow >= 2 && twopow % 2 == 0) {
			fftt.N1 = 4;
		} else {
			int threepow = 0;
			int n = N;
			while (n % 3 == 0) {
				threepow++;
				n /= 3;
			}
			if (false && threepow >= 2 && threepow % 2 == 0) {
				fftt.N1 = 9;
			} else {
				bool found = false;
				for (int r = 2; r <= rader_radix; r++) {
					if (is_prime(r) && N % r == 0) {
						fftt.N1 = r;
						found = true;
						break;
					}
				}
				if (!found) {
					fftt.type = RADERS;
					fftt.N1 = N;
				}
			}
		}
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
