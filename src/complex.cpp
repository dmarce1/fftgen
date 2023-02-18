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

void print_complex_short_fft(int r, std::vector<std::string> in_, std::vector<std::string> out_) {
	std::vector<const char*> in(in_.size());
	std::vector<const char*> out(out_.size());
	for (int i = 0; i < in.size(); i++) {
		in[i] = in_[i].c_str();
		out[i] = out_[i].c_str();
	}
	switch (r) {
	case 2:
		print("const auto tr = %s;\n", in[0]);
		print("const auto ti = %s;\n", in[1]);
		print("%s = tr + %s;\n", out[0], in[2]);
		print("%s = ti + %s;\n", out[1], in[3]);
		print("%s = tr - %s;\n", out[2], in[2]);
		print("%s = ti - %s;\n", out[3], in[3]);
		break;
	case 3:
		print("const auto tr1 = %s + %s;\n", in[2], in[4]);
		print("const auto ti1 = %s + %s;\n", in[3], in[5]);
		print("const auto tr2 = std::fma(tr1, -0.5, %s);\n", in[0]);
		print("const auto ti2 = std::fma(ti1, -0.5, %s);\n", in[1]);
		print("const auto tr3 = (%24.17e) * (%s - %s);\n", sqrt(3) * 0.5, in[2], in[4]);
		print("const auto ti3 = (%24.17e) * (%s - %s);\n", sqrt(3) * 0.5, in[3], in[5]);
		print("%s = %s + tr1;\n", out[0], in[0]);
		print("%s = %s + ti1;\n", out[1], in[1]);
		print("%s = tr2 + ti3;\n", out[2]);
		print("%s = ti2 - tr3;\n", out[3]);
		print("%s = tr2 - ti3;\n", out[4]);
		print("%s = ti2 + tr3;\n", out[5]);
		break;
	case 4:
		print("const auto tr1 = %s + %s;\n", in[0], in[4]);
		print("const auto ti1 = %s + %s;\n", in[1], in[5]);
		print("const auto tr3 = %s - %s;\n", in[0], in[4]);
		print("const auto ti3 = %s - %s;\n", in[1], in[5]);
		print("const auto tr2 = %s + %s;\n", in[2], in[6]);
		print("const auto ti2 = %s + %s;\n", in[3], in[7]);
		print("const auto tr4 = %s - %s;\n", in[2], in[6]);
		print("const auto ti4 = %s - %s;\n", in[3], in[7]);
		print("%s = tr1 + tr2;\n", out[0]);
		print("%s = ti1 + ti2;\n", out[1]);
		print("%s = tr3 + ti4;\n", out[2]);
		print("%s = ti3 - tr4;\n", out[3]);
		print("%s = tr1 - tr2;\n", out[4]);
		print("%s = ti1 - ti2;\n", out[5]);
		print("%s = tr3 - ti4;\n", out[6]);
		print("%s = ti3 + tr4;\n", out[7]);
		break;
	case 5:
		print("const auto tr1 = %s + %s;\n", in[2], in[8]);
		print("const auto ti1 = %s + %s;\n", in[3], in[9]);
		print("const auto tr2 = %s + %s;\n", in[4], in[6]);
		print("const auto ti2 = %s + %s;\n", in[5], in[7]);
		print("const auto tr3 = %s - %s;\n", in[2], in[8]);
		print("const auto ti3 = %s - %s;\n", in[3], in[9]);
		print("const auto tr4 = %s - %s;\n", in[4], in[6]);
		print("const auto ti4 = %s - %s;\n", in[5], in[7]);
		print("const auto tr5 = tr1 + tr2;\n");
		print("const auto ti5 = ti1 + ti2;\n");
		print("const auto tr6 = (%24.17e) * (tr1 - tr2);\n", sqrt(5) * 0.25);
		print("const auto ti6 = (%24.17e) * (ti1 - ti2);\n", sqrt(5) * 0.25);
		print("const auto tr7 = std::fma(tr5, -0.25, %s);\n", in[0]);
		print("const auto ti7 = std::fma(ti5, -0.25, %s);\n", in[1]);
		print("const auto tr8 = tr7 + tr6;\n");
		print("const auto ti8 = ti7 + ti6;\n");
		print("const auto tr9 = tr7 - tr6;\n");
		print("const auto ti9 = ti7 - ti6;\n");
		print("const auto tr10 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
		print("const auto ti10 = std::fma((%24.17e), ti3, (%24.17e) * ti4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
		print("const auto tr11 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
		print("const auto ti11 = std::fma((%24.17e), ti3, (%24.17e) * ti4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
		print("%s = %s + tr5;\n", out[0], in[0]);
		print("%s = %s + ti5;\n", out[1], in[1]);
		print("%s = tr8 + ti10;\n", out[2]);
		print("%s = ti8 - tr10;\n", out[3]);
		print("%s = tr9 + ti11;\n", out[4]);
		print("%s = ti9 - tr11;\n", out[5]);
		print("%s = tr9 - ti11;\n", out[6]);
		print("%s = ti9 + tr11;\n", out[7]);
		print("%s = tr8 - ti10;\n", out[8]);
		print("%s = ti8 + tr10;\n", out[9]);
		break;
	case 6:
		print("const auto tr1 = %s + %s;\n", in[4], in[8]);
		print("const auto ti1 = %s + %s;\n", in[5], in[9]);
		print("const auto tr2 = std::fma(tr1, -0.5, %s);\n", in[0]);
		print("const auto ti2 = std::fma(ti1, -0.5, %s);\n", in[1]);
		print("const auto tr3 = (%24.17e) * (%s - %s);\n", sin(M_PI / 3.0), in[4], in[8]);
		print("const auto ti3 = (%24.17e) * (%s - %s);\n", sin(M_PI / 3.0), in[5], in[9]);
		print("const auto tr4 = %s + %s;\n", in[10], in[2]);
		print("const auto ti4 = %s + %s;\n", in[11], in[3]);
		print("const auto tr5 = std::fma(tr4, -0.5, %s);\n", in[6]);
		print("const auto ti5 = std::fma(ti4, -0.5, %s);\n", in[7]);
		print("const auto tr6 = (%24.17e) * (%s - %s);\n", sin(M_PI / 3.0), in[10], in[2]);
		print("const auto ti6 = (%24.17e) * (%s - %s);\n", sin(M_PI / 3.0), in[11], in[3]);
		print("const auto tr7 = %s + tr1;\n", in[0]);
		print("const auto ti7 = %s + ti1;\n", in[1]);
		print("const auto tr8 = tr2 + ti3;\n");
		print("const auto ti8 = ti2 - tr3;\n");
		print("const auto tr9 = tr2 - ti3;\n");
		print("const auto ti9 = ti2 + tr3;\n");
		print("const auto tr10 = %s + tr4;\n", in[6]);
		print("const auto ti10 = %s + ti4;\n", in[7]);
		print("const auto tr11 = tr5 + ti6;\n");
		print("const auto ti11 = ti5 - tr6;\n");
		print("const auto tr12 = tr5 - ti6;\n");
		print("const auto ti12 = ti5 + tr6;\n");
		print("%s = tr7 + tr10;\n", out[0]);
		print("%s = ti7 + ti10;\n", out[1]);
		print("%s = tr8 - tr11;\n", out[2]);
		print("%s = ti8 - ti11;\n", out[3]);
		print("%s = tr9 + tr12;\n", out[4]);
		print("%s = ti9 + ti12;\n", out[5]);
		print("%s = tr7 - tr10;\n", out[6]);
		print("%s = ti7 - ti10;\n", out[7]);
		print("%s = ti8 + ti11;\n", out[9]);
		print("%s = tr8 + tr11;\n", out[8]);
		print("%s = tr9 - tr12;\n", out[10]);
		print("%s = ti9 - ti12;\n", out[11]);
		break;
	case 7: {
		constexpr double theta = 2.0 * M_PI / 7.0;
		const double c1 = cos(theta);
		const double c2 = cos(2.0 * theta);
		const double c3 = cos(3.0 * theta);
		const double c4 = sin(theta);
		const double c5 = sin(2.0 * theta);
		const double c6 = sin(3.0 * theta);
		print("const auto tr1 = %s + %s;\n", in[2], in[12]);
		print("const auto ti1 = %s + %s;\n", in[3], in[13]);
		print("const auto tr2 = %s + %s;\n", in[4], in[10]);
		print("const auto ti2 = %s + %s;\n", in[5], in[11]);
		print("const auto tr3 = %s + %s;\n", in[6], in[8]);
		print("const auto ti3 = %s + %s;\n", in[7], in[9]);
		print("const auto tr4 = %s - %s;\n", in[2], in[12]);
		print("const auto ti4 = %s - %s;\n", in[3], in[13]);
		print("const auto tr5 = %s - %s;\n", in[4], in[10]);
		print("const auto ti5 = %s - %s;\n", in[5], in[11]);
		print("const auto tr6 = %s - %s;\n", in[6], in[8]);
		print("const auto ti6 = %s - %s;\n", in[7], in[9]);
		print("const auto tr7 = %s - 0.5 * tr3;\n", in[0]);
		print("const auto ti7 = %s - 0.5 * ti3;\n", in[1]);
		print("const auto tr8 = tr1 - tr3;\n");
		print("const auto ti8 = ti1 - ti3;\n");
		print("const auto tr9 = tr2 - tr3;\n");
		print("const auto ti9 = ti2 - ti3;\n");
		print("const auto yr0 = %s + tr1 + tr2 + tr3;\n", in[0]);
		print("const auto yi0 = %s + ti1 + ti2 + ti3;\n", in[1]);
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
		print("%s = yr0;\n", out[0]);
		print("%s = yi0;\n", out[1]);
		print("%s = yr1 + yi6;\n", out[2]);
		print("%s = yi1 - yr6;\n", out[3]);
		print("%s = yr2 + yi5;\n", out[4]);
		print("%s = yi2 - yr5;\n", out[5]);
		print("%s = yr3 + yi4;\n", out[6]);
		print("%s = yi3 - yr4;\n", out[7]);
		print("%s = yr3 - yi4;\n", out[8]);
		print("%s = yi3 + yr4;\n", out[9]);
		print("%s = yr2 - yi5;\n", out[10]);
		print("%s = yi2 + yr5;\n", out[11]);
		print("%s = yr1 - yi6;\n", out[12]);
		print("%s = yi1 + yr6;\n", out[13]);
	}
		break;
	case 8:
		print("const auto tr1 = %s + %s;\n", in[0], in[8]);
		print("const auto ti1 = %s + %s;\n", in[1], in[9]);
		print("const auto tr2 = %s - %s;\n", in[0], in[8]);
		print("const auto ti2 = %s - %s;\n", in[1], in[9]);
		print("const auto tr3 = %s + %s;\n", in[2], in[10]);
		print("const auto ti3 = %s + %s;\n", in[3], in[11]);
		print("const auto tr4 = %s - %s;\n", in[2], in[10]);
		print("const auto ti4 = %s - %s;\n", in[3], in[11]);
		print("const auto tr5 = %s + %s;\n", in[4], in[12]);
		print("const auto ti5 = %s + %s;\n", in[5], in[13]);
		print("const auto tr6 = %s - %s;\n", in[4], in[12]);
		print("const auto ti6 = %s - %s;\n", in[5], in[13]);
		print("const auto tr7 = %s + %s;\n", in[6], in[14]);
		print("const auto ti7 = %s + %s;\n", in[7], in[15]);
		print("const auto tr8 = %s - %s;\n", in[6], in[14]);
		print("const auto ti8 = %s - %s;\n", in[7], in[15]);
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
		print("%s = yr0;\n", out[0]);
		print("%s = yi0;\n", out[1]);
		print("%s = yr1 + yi7;\n", out[2]);
		print("%s = yi1 - yr7;\n", out[3]);
		print("%s = yr2 + yi6;\n", out[4]);
		print("%s = yi2 - yr6;\n", out[5]);
		print("%s = yr3 + yi5;\n", out[6]);
		print("%s = yi3 - yr5;\n", out[7]);
		print("%s = yr4;\n", out[8]);
		print("%s = yi4;\n", out[9]);
		print("%s = yr3 - yi5;\n", out[10]);
		print("%s = yi3 + yr5;\n", out[11]);
		print("%s = yr2 - yi6;\n", out[12]);
		print("%s = yi2 + yr6;\n", out[13]);
		print("%s = yr1 - yi7;\n", out[14]);
		print("%s = yi1 + yr7;\n", out[15]);
		break;
	case 9: {
		const double u = -2.0 * M_PI / 9.0;
		const double ar0 = 1.0;
		const double ar1 = -0.5;
		const double ai2 = sin(3.0 * u);
		const double ar3 = cos(3.0 * u) - 1.0;
		const double ai4 = sin(3.0 * u);
		const double ar5 = (2.0 * cos(u) - cos(2.0 * u) - cos(4.0 * u)) / 3.0;
		const double ar6 = (cos(u) + cos(2.0 * u) - 2.0 * cos(4.0 * u)) / 3.0;
		const double ar7 = (cos(u) - 2.0 * cos(2.0 * u) + cos(4.0 * u)) / 3.0;
		const double ai8 = (2.0 * sin(u) + sin(2.0 * u) - sin(4.0 * u)) / 3.0;
		const double ai9 = (sin(u) - sin(2.0 * u) - 2.0 * sin(4.0 * u)) / 3.0;
		const double ai10 = (sin(u) + 2.0 * sin(2.0 * u) + sin(4.0 * u)) / 3.0;
		print("const auto sr1 = %s + %s;\n", in[2], in[16]);
		print("const auto si1 = %s + %s;\n", in[3], in[17]);
		print("const auto sr2 = %s - %s;\n", in[2], in[16]);
		print("const auto si2 = %s - %s;\n", in[3], in[17]);
		print("const auto sr3 = %s + %s;\n", in[14], in[4]);
		print("const auto si3 = %s + %s;\n", in[15], in[5]);
		print("const auto sr4 = %s - %s;\n", in[14], in[4]);
		print("const auto si4 = %s - %s;\n", in[15], in[5]);
		print("const auto sr5 = %s + %s;\n", in[6], in[12]);
		print("const auto si5 = %s + %s;\n", in[7], in[13]);
		print("const auto sr6 = %s - %s;\n", in[6], in[12]);
		print("const auto si6 = %s - %s;\n", in[7], in[13]);
		print("const auto sr7 = %s + %s;\n", in[8], in[10]);
		print("const auto si7 = %s + %s;\n", in[9], in[11]);
		print("const auto sr8 = %s - %s;\n", in[8], in[10]);
		print("const auto si8 = %s - %s;\n", in[9], in[11]);
		print("const auto sr9 = sr1 + sr3;\n");
		print("const auto si9 = si1 + si3;\n");
		print("const auto sr10 = sr9 + sr7;\n");
		print("const auto si10 = si9 + si7;\n");
		print("const auto sr11 = sr10 + sr5;\n");
		print("const auto si11 = si10 + si5;\n");
		print("const auto sr12 = sr11 + %s;\n", in[0]);
		print("const auto si12 = si11 + %s;\n", in[1]);
		print("const auto sr13 = sr2 + sr4;\n");
		print("const auto si13 = si2 + si4;\n");
		print("const auto sr14 = sr13 + sr8;\n");
		print("const auto si14 = si13 + si8;\n");
		print("const auto sr15 = sr1 - sr3;\n");
		print("const auto si15 = si1 - si3;\n");
		print("const auto sr16 = sr3 - sr7;\n");
		print("const auto si16 = si3 - si7;\n");
		print("const auto sr17 = sr7 - sr1;\n");
		print("const auto si17 = si7 - si1;\n");
		print("const auto sr18 = sr2 - sr4;\n");
		print("const auto si18 = si2 - si4;\n");
		print("const auto sr19 = sr4 - sr8;\n");
		print("const auto si19 = si4 - si8;\n");
		print("const auto sr20 = sr8 - sr2;\n");
		print("const auto si20 = si8 - si2;\n");
		print("const auto mr0 = sr12;\n");
		print("const auto mi0 = si12;\n");
		print("const auto mr1 = sr10 * (%.17e);\n", ar1);
		print("const auto mi1 = si10 * (%.17e);\n", ar1);
		print("const auto mr2 = si14 * (%.17e);\n", -ai2);
		print("const auto mi2 = sr14 * (%.17e);\n", ai2);
		print("const auto mr3 = sr5 * (%.17e);\n", ar3);
		print("const auto mi3 = si5 * (%.17e);\n", ar3);
		print("const auto mr4 = si6 * (%.17e);\n", -ai4);
		print("const auto mi4 = sr6 * (%.17e);\n", ai4);
		print("const auto mr5 = sr15 * (%.17e);\n", ar5);
		print("const auto mi5 = si15 * (%.17e);\n", ar5);
		print("const auto mr6 = sr16 * (%.17e);\n", ar6);
		print("const auto mi6 = si16 * (%.17e);\n", ar6);
		print("const auto mr7 = sr17 * (%.17e);\n", ar7);
		print("const auto mi7 = si17 * (%.17e);\n", ar7);
		print("const auto mr8 = si18 * (%.17e);\n", -ai8);
		print("const auto mi8 = sr18 * (%.17e);\n", ai8);
		print("const auto mr9 = si19 * (%.17e);\n", -ai9);
		print("const auto mi9 = sr19 * (%.17e);\n", ai9);
		print("const auto mr10 = si20 * (%.17e);\n", -ai10);
		print("const auto mi10 = sr20 * (%.17e);\n", ai10);
		print("const auto sr21 = mr1 + mr1;\n");
		print("const auto si21 = mi1 + mi1;\n");
		print("const auto sr22 = sr21 + mr1;\n");
		print("const auto si22 = si21 + mi1;\n");
		print("const auto sr23 = sr22 + mr0;\n");
		print("const auto si23 = si22 + mi0;\n");
		print("const auto sr24 = sr23 + mr2;\n");
		print("const auto si24 = si23 + mi2;\n");
		print("const auto sr25 = sr23 - mr2;\n");
		print("const auto si25 = si23 - mi2;\n");
		print("const auto sr26 = mr0 + mr3;\n");
		print("const auto si26 = mi0 + mi3;\n");
		print("const auto sr27 = sr26 + sr21;\n");
		print("const auto si27 = si26 + si21;\n");
		print("const auto sr28 = sr27 + mr5;\n");
		print("const auto si28 = si27 + mi5;\n");
		print("const auto sr29 = sr28 + mr6;\n");
		print("const auto si29 = si28 + mi6;\n");
		print("const auto sr30 = sr27 - mr6;\n");
		print("const auto si30 = si27 - mi6;\n");
		print("const auto sr31 = sr30 + mr7;\n");
		print("const auto si31 = si30 + mi7;\n");
		print("const auto sr32 = sr27 - mr5;\n");
		print("const auto si32 = si27 - mi5;\n");
		print("const auto sr33 = sr32 - mr7;\n");
		print("const auto si33 = si32 - mi7;\n");
		print("const auto sr34 = mr4 + mr8;\n");
		print("const auto si34 = mi4 + mi8;\n");
		print("const auto sr35 = sr34 + mr9;\n");
		print("const auto si35 = si34 + mi9;\n");
		print("const auto sr36 = mr4 - mr9;\n");
		print("const auto si36 = mi4 - mi9;\n");
		print("const auto sr37 = sr36 + mr10;\n");
		print("const auto si37 = si36 + mi10;\n");
		print("const auto sr38 = mr4 - mr8;\n");
		print("const auto si38 = mi4 - mi8;\n");
		print("const auto sr39 = sr38 - mr10;\n");
		print("const auto si39 = si38 - mi10;\n");
		print("const auto sr40 = sr29 + sr35;\n");
		print("const auto si40 = si29 + si35;\n");
		print("const auto sr41 = sr29 - sr35;\n");
		print("const auto si41 = si29 - si35;\n");
		print("const auto sr42 = sr31 + sr37;\n");
		print("const auto si42 = si31 + si37;\n");
		print("const auto sr43 = sr31 - sr37;\n");
		print("const auto si43 = si31 - si37;\n");
		print("const auto sr44 = sr33 + sr39;\n");
		print("const auto si44 = si33 + si39;\n");
		print("const auto sr45 = sr33 - sr39;\n");
		print("const auto si45 = si33 - si39;\n");
		print("%s = mr0;\n", out[0]);
		print("%s = mi0;\n", out[1]);
		print("%s = sr40;\n", out[2]);
		print("%s = si40;\n", out[3]);
		print("%s = sr43;\n", out[4]);
		print("%s = si43;\n", out[5]);
		print("%s = sr24;\n", out[6]);
		print("%s = si24;\n", out[7]);
		print("%s = sr44;\n", out[8]);
		print("%s = si44;\n", out[9]);
		print("%s = sr45;\n", out[10]);
		print("%s = si45;\n", out[11]);
		print("%s = sr25;\n", out[12]);
		print("%s = si25;\n", out[13]);
		print("%s = sr42;\n", out[14]);
		print("%s = si42;\n", out[15]);
		print("%s = sr41;\n", out[16]);
		print("%s = si41;\n", out[17]);

	}
		break;
	case 16: {
		const auto theta = M_PI / 8.0;
		const auto c2 = cos(theta);
		const auto c3 = sin(theta);
		const auto c4 = cos(2.0 * theta);
		const auto c5 = c4;
		const auto c6 = c3;
		const auto c7 = c2;
		print("const auto tr1 = %s + %s;\n", in[2 * 0 + 0], in[2 * 8 + 0]);
		print("const auto ti1 = %s + %s;\n", in[2 * 0 + 1], in[2 * 8 + 1]);
		print("const auto tr2 = %s + %s;\n", in[2 * 4 + 0], in[2 * 12 + 0]);
		print("const auto ti2 = %s + %s;\n", in[2 * 4 + 1], in[2 * 12 + 1]);
		print("const auto tr3 = %s - %s;\n", in[2 * 0 + 0], in[2 * 8 + 0]);
		print("const auto ti3 = %s - %s;\n", in[2 * 0 + 1], in[2 * 8 + 1]);
		print("const auto tr4 = %s - %s;\n", in[2 * 4 + 0], in[2 * 12 + 0]);
		print("const auto ti4 = %s - %s;\n", in[2 * 4 + 1], in[2 * 12 + 1]);
		print("const auto tr5 = tr1 + tr2;\n");
		print("const auto ti5 = ti1 + ti2;\n");
		print("const auto tr6 = tr1 - tr2;\n");
		print("const auto ti6 = ti1 - ti2;\n");
		print("const auto tr7 = %s + %s;\n", in[2 * 1 + 0], in[2 * 9 + 0]);
		print("const auto ti7 = %s + %s;\n", in[2 * 1 + 1], in[2 * 9 + 1]);
		print("const auto tr8 = %s + %s;\n", in[2 * 5 + 0], in[2 * 13 + 0]);
		print("const auto ti8 = %s + %s;\n", in[2 * 5 + 1], in[2 * 13 + 1]);
		print("const auto tr9 = %s - %s;\n", in[2 * 1 + 0], in[2 * 9 + 0]);
		print("const auto ti9 = %s - %s;\n", in[2 * 1 + 1], in[2 * 9 + 1]);
		print("const auto tr10 = %s - %s;\n", in[2 * 5 + 0], in[2 * 13 + 0]);
		print("const auto ti10 = %s - %s;\n", in[2 * 5 + 1], in[2 * 13 + 1]);
		print("const auto tr11 = tr7 + tr8;\n");
		print("const auto ti11 = ti7 + ti8;\n");
		print("const auto tr12 = tr7 - tr8;\n");
		print("const auto ti12 = ti7 - ti8;\n");
		print("const auto tr13 = %s + %s;\n", in[2 * 2 + 0], in[2 * 10 + 0]);
		print("const auto ti13 = %s + %s;\n", in[2 * 2 + 1], in[2 * 10 + 1]);
		print("const auto tr14 = %s + %s;\n", in[2 * 6 + 0], in[2 * 14 + 0]);
		print("const auto ti14 = %s + %s;\n", in[2 * 6 + 1], in[2 * 14 + 1]);
		print("const auto tr15 = %s - %s;\n", in[2 * 2 + 0], in[2 * 10 + 0]);
		print("const auto ti15 = %s - %s;\n", in[2 * 2 + 1], in[2 * 10 + 1]);
		print("const auto tr16 = %s - %s;\n", in[2 * 6 + 0], in[2 * 14 + 0]);
		print("const auto ti16 = %s - %s;\n", in[2 * 6 + 1], in[2 * 14 + 1]);
		print("const auto tr17 = tr13 + tr14;\n");
		print("const auto ti17 = ti13 + ti14;\n");
		print("const auto tr18 = (%.17e) * (tr15 - tr16);\n", c4);
		print("const auto ti18 = (%.17e) * (ti15 - ti16);\n", c4);
		print("const auto tr19 = (%.17e) * (tr15 + tr16);\n", c5);
		print("const auto ti19 = (%.17e) * (ti15 + ti16);\n", c5);
		print("const auto tr20 = tr13 - tr14;\n");
		print("const auto ti20 = ti13 - ti14;\n");
		print("const auto tr21 = %s + %s;\n", in[2 * 3 + 0], in[2 * 11 + 0]);
		print("const auto ti21 = %s + %s;\n", in[2 * 3 + 1], in[2 * 11 + 1]);
		print("const auto tr22 = %s + %s;\n", in[2 * 7 + 0], in[2 * 15 + 0]);
		print("const auto ti22 = %s + %s;\n", in[2 * 7 + 1], in[2 * 15 + 1]);
		print("const auto tr23 = %s - %s;\n", in[2 * 3 + 0], in[2 * 11 + 0]);
		print("const auto ti23 = %s - %s;\n", in[2 * 3 + 1], in[2 * 11 + 1]);
		print("const auto tr24 = %s - %s;\n", in[2 * 7 + 0], in[2 * 15 + 0]);
		print("const auto ti24 = %s - %s;\n", in[2 * 7 + 1], in[2 * 15 + 1]);
		print("const auto tr25 = tr21 + tr22;\n");
		print("const auto ti25 = ti21 + ti22;\n");
		print("const auto tr26 = tr21 - tr22;\n");
		print("const auto ti26 = ti21 - ti22;\n");
		print("const auto tr27 = tr9 + tr24;\n");
		print("const auto ti27 = ti9 + ti24;\n");
		print("const auto tr28 = tr10 + tr23;\n");
		print("const auto ti28 = ti10 + ti23;\n");
		print("const auto tr29 = tr9 - tr24;\n");
		print("const auto ti29 = ti9 - ti24;\n");
		print("const auto tr30 = tr10 - tr23;\n");
		print("const auto ti30 = ti10 - ti23;\n");
		print("const auto tr31 = tr5 + tr17;\n");
		print("const auto ti31 = ti5 + ti17;\n");
		print("const auto tr32 = tr11 + tr25;\n");
		print("const auto ti32 = ti11 + ti25;\n");
		print("const auto tr33 = tr3 + tr18;\n");
		print("const auto ti33 = ti3 + ti18;\n");
		print("const auto tr34 = (%.17e) * tr29 + (%.17e) * tr30;\n", c2, -c6);
		print("const auto ti34 = (%.17e) * ti29 + (%.17e) * ti30;\n", c2, -c6);
		print("const auto tr35 = tr3 - tr18;\n");
		print("const auto ti35 = ti3 - ti18;\n");
		print("const auto tr36 = (%.17e) * tr27 + (%.17e) * tr28;\n", c7, -c3);
		print("const auto ti36 = (%.17e) * ti27 + (%.17e) * ti28;\n", c7, -c3);
		print("const auto tr37 = tr4 + tr19;\n");
		print("const auto ti37 = ti4 + ti19;\n");
		print("const auto tr38 = (%.17e) * tr27 + (%.17e) * tr28;\n", c3, c7);
		print("const auto ti38 = (%.17e) * ti27 + (%.17e) * ti28;\n", c3, c7);
		print("const auto tr39 = tr4 - tr19;\n");
		print("const auto ti39 = ti4 - ti19;\n");
		print("const auto tr40 = (%.17e) * tr29 + (%.17e) * tr30;\n", c6, c2);
		print("const auto ti40 = (%.17e) * ti29 + (%.17e) * ti30;\n", c6, c2);
		print("const auto tr41 = (%.17e) * (tr12 - tr26);\n", c4);
		print("const auto ti41 = (%.17e) * (ti12 - ti26);\n", c4);
		print("const auto tr42 = (%.17e) * (tr12 + tr26);\n", c5);
		print("const auto ti42 = (%.17e) * (ti12 + ti26);\n", c5);
		print("const auto yr0 = tr31 + tr32;\n");
		print("const auto yi0 = ti31 + ti32;\n");
		print("const auto yr1 = tr33 + tr34;\n");
		print("const auto yi1 = ti33 + ti34;\n");
		print("const auto yr2 = tr6 + tr41;\n");
		print("const auto yi2 = ti6 + ti41;\n");
		print("const auto yr3 = tr35 + tr40;\n");
		print("const auto yi3 = ti35 + ti40;\n");
		print("const auto yr4 = tr5 - tr17;\n");
		print("const auto yi4 = ti5 - ti17;\n");
		print("const auto yr5 = tr35 - tr40;\n");
		print("const auto yi5 = ti35 - ti40;\n");
		print("const auto yr6 = tr6 - tr41;\n");
		print("const auto yi6 = ti6 - ti41;\n");
		print("const auto yr7 = tr33 - tr34;\n");
		print("const auto yi7 = ti33 - ti34;\n");
		print("const auto yr8 = tr31 - tr32;\n");
		print("const auto yi8 = ti31 - ti32;\n");
		print("const auto yr9 = tr38 - tr37;\n");
		print("const auto yi9 = ti38 - ti37;\n");
		print("const auto yr10 = tr42 - tr20;\n");
		print("const auto yi10 = ti42 - ti20;\n");
		print("const auto yr11 = tr36 + tr39;\n");
		print("const auto yi11 = ti36 + ti39;\n");
		print("const auto yr12 = tr11 - tr25;\n");
		print("const auto yi12 = ti11 - ti25;\n");
		print("const auto yr13 = tr36 - tr39;\n");
		print("const auto yi13 = ti36 - ti39;\n");
		print("const auto yr14 = tr42 + tr20;\n");
		print("const auto yi14 = ti42 + ti20;\n");
		print("const auto yr15 = tr38 + tr37;\n");
		print("const auto yi15 = ti38 + ti37;\n");
		print("%s = yr0;\n", out[2 * 0 + 0]);
		print("%s = yi0;\n", out[2 * 0 + 1]);
		for (int i = 1; i <= 7; i++) {
			print("%s = yr%i + yi%i;\n", out[2 * i + 0], i, 16 - i);
			print("%s = yi%i - yr%i;\n", out[2 * i + 1], i, 16 - i);
		}
		print("%s = yr8;\n", out[2 * 8 + 0]);
		print("%s = yi8;\n", out[2 * 8 + 1]);
		for (int i = 1; i <= 7; i++) {
			print("%s = yr%i - yi%i;\n", out[2 * (16 - i) + 0], i, 16 - i);
			print("%s = yi%i + yr%i;\n", out[2 * (16 - i) + 1], i, 16 - i);
		}
	}
		break;
	case 10:
	case 12:
	case 14:
	case 15: {
		print("double zr0, ");
		for (int n = 1; n < r - 1; n++) {
			print_notab("zr%i, ", n);
		}
		print_notab("zr%i;\n", r - 1);
		print("double zi0, ");
		for (int n = 1; n < r - 1; n++) {
			print_notab("zi%i, ", n);
		}
		print_notab("zi%i;\n", r - 1);
		int N1;
		int N2;
		if (r == 10) {
			N1 = 2;
			N2 = 5;
		} else if (r == 12) {
			N1 = 3;
			N2 = 4;
		} else if (r == 14) {
			N1 = 2;
			N2 = 7;
		} else if (r == 15) {
			N1 = 3;
			N2 = 5;
		}
		std::vector<std::vector<std::string>> out1(N2, std::vector<std::string>(2 * N1));
		std::vector<std::vector<std::string>> in1(N2, std::vector<std::string>(2 * N1));
		std::vector<std::vector<std::string>> out2(N1, std::vector<std::string>(2 * N2));
		std::vector<std::vector<std::string>> in2(N1, std::vector<std::string>(2 * N2));
		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				const int nn = ((n1 * N2 + n2 * N1) % r);
				in1[n2][2 * n1] = in[2 * nn];
				in1[n2][2 * n1 + 1] = in[2 * nn + 1];
				out1[n2][2 * n1] = std::string("zr") + std::to_string(n1 + n2 * N1);
				out1[n2][2 * n1 + 1] = std::string("zi") + std::to_string(n1 + n2 * N1);
				in2[n1][2 * n2] = std::string("zr") + std::to_string(n1 + n2 * N1);
				in2[n1][2 * n2 + 1] = std::string("zi") + std::to_string(n1 + n2 * N1);
			}
		}
		for (int n = 0; n < r; n++) {
			int n1 = n % N1;
			int n2 = n % N2;
			out2[n1][2 * n2] = out[2 * n];
			out2[n1][2 * n2 + 1] = out[2 * n + 1];
		}
		for (int n2 = 0; n2 < N2; n2++) {
			print("{\n");
			indent();
			print_complex_short_fft(N1, in1[n2], out1[n2]);
			deindent();
			print("}\n");
		}
		for (int n1 = 0; n1 < N1; n1++) {
			print("{\n");
			indent();
			print_complex_short_fft(N2, in2[n1], out2[n1]);
			deindent();
			print("}\n");
		}
	}

		break;
	default:
		for (int j = 1; j <= (r - 1) / 2; j++) {
			print("const auto txp%i = %s + %s;\n", j, in[2 * j], in[2 * (r - j)]);
		}
		for (int j = 1; j <= (r - 1) / 2; j++) {
			print("const auto txm%i = %s - %s;\n", j, in[2 * j], in[2 * (r - j)]);
		}
		for (int j = 1; j <= (r - 1) / 2; j++) {
			print("const auto typ%i = %s + %s;\n", j, in[2 * j + 1], in[2 * (r - j) + 1]);
		}
		for (int j = 1; j <= (r - 1) / 2; j++) {
			print("const auto tym%i = %s - %s;\n", j, in[2 * j + 1], in[2 * (r - j) + 1]);
		}
		for (int i = 1; i <= (r - 1) / 2; i++) {
			print("auto ap%i = %s;\n", i, in[0]);
			print("auto bp%i = %s;\n", i, in[1]);
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
		print("%s = ", out[0]);
		for (int i = 0; i < r; i++) {
			print_notab(" + %s", in[2 * i]);
		}
		print_notab(";\n");
		print("%s = ", out[1]);
		for (int i = 0; i < r; i++) {
			print_notab(" + %s", in[2 * i + 1]);
		}
		print_notab(";\n");
		for (int i = 1; i <= (r - 1) / 2; i++) {
			print("%s = ap%i + am%i;\n", out[2 * i], i, i);
			print("%s = bp%i - bm%i;\n", out[2 * i + 1], i, i);
			print("%s = ap%i - am%i;\n", out[2 * (r - i)], i, i);
			print("%s = bp%i + bm%i;\n", out[2 * (r - i) + 1], i, i);
		}
	}

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

fft_type best_radix(int N, int o, bool first) {
	fft_type fftt;
	fftt.type = RADIX;
	fftt.nops = 0;
	const int rader_radix = is_prime(N) ? RADER_LEN : 1000000000;
	if (N % 14 == 0) {
		fftt.N1 = 14;
	} else if (N % 9 == 0) {
		fftt.N1 = 9;
	} else if (N % 15 == 0) {
		fftt.N1 = 15;
	} else if (N % 10 == 0) {
		fftt.N1 = 10;
	} else if (N % 12 == 0) {
		fftt.N1 = 12;
	} else {
		if (N % 6 == 0) {
			fftt.N1 = 6;
		} else {
			int twopow = 0;
			int n = N;
			while (n % 2 == 0) {
				twopow++;
				n /= 2;
			}
			if (false && twopow >= 4 && twopow % 4 == 0) {
				fftt.N1 = 16;
			} else if (twopow >= 3 && twopow % 3 == 0) {
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
