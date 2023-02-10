#include "fftgen.hpp"
#include <list>

void print_z(int zi, int twi, int loci, int k, int r, int N, int o, bool ref) {
	const auto W = twiddle(twi * k, N);
	const int i = k + loci * N / r;
	if (twi * k == 0) {
		print("const auto%s zr%i = x[%i];\n", ref ? "&" : "", zi, index(o, i, 0, N));
		print("const auto%s zi%i = x[%i];\n", ref ? "&" : "", zi, index(o, i, 1, N));
	} else if (twi * k == N / 8 && N % 8 == 0) {
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
	}
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
	case 4: {
		std::vector<int> J;
		for (int k = 0; k < N / 2; k++) {
			J.push_back(I[2 * k]);
		}
		auto K = fft_bitr(N / 2, o, J);
		L.insert(L.end(), K.begin(), K.end());
	}
		for (int n = 1; n < 4; n += 2) {
			std::vector<int> J;
			for (int k = 0; k < N / r; k++) {
				J.push_back(I[n + 4 * k]);
			}
			auto K = fft_bitr(N / 4, o + (2 + (n / 2)) * N / 4, J);
			L.insert(L.end(), K.begin(), K.end());
		}
		break;
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
	} else {
		return indices;
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
	if (r == 4) {
		fft(N / 2, o);
		fft(N / 4, o + N / 2);
		fft(N / 4, o + 3 * N / 4);
	} else {
		for (int n = 0; n < r; n++) {
			fft(N / r, o + n * N / r);
		}
	}

	for (int k = 0; k < N / r; k++) {
		print("{\n");
		indent();
		if (r != 4) {
			print_z(0, 0, 0, k, r, N, o, r != 2);
			print_z(1, 1, 1, k, r, N, o, true);
			for (int i = 2; i < r; i++) {
				print_z(i, i, i, k, r, N, o, r != 2);
			}
		} else {
			print_z(0, 0, 0, k, r, N, o, false);
			print_z(1, 1, 2, k, r, N, o, true);
			print_z(2, 0, 1, k, r, N, o, false);
			print_z(3, 3, 3, k, r, N, o, true);
		}

		switch (r) {
		case 1:
			return;
		case 2:
			print("x[%i] = zr0 + zr1;\n", index(o, k, 0, N));
			print("x[%i] = zi0 + zi1;\n", index(o, k, 1, N));
			print("x[%i] = zr0 - zr1;\n", index(o, k + N / 2, 0, N));
			print("x[%i] = zi0 - zi1;\n", index(o, k + N / 2, 1, N));
			break;
		case 3:
			print("const auto tr1 = zr1 + zr2;\n");
			print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
			print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
			print("const auto ti1 = zi1 + zi2;\n");
			print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
			print("const auto ti3 = (%24.17e) * (zi1 - zi2);\n", sqrt(3) * 0.5);
			print("x[%i] = zr0 + tr1;\n", index(o, k, 0, N));
			print("x[%i] = zi0 + ti1;\n", index(o, k, 1, N));
			print("x[%i] = tr2 + ti3;\n", index(o, k + N / 3, 0, N));
			print("x[%i] = ti2 - tr3;\n", index(o, k + N / 3, 1, N));
			print("x[%i] = tr2 - ti3;\n", index(o, k + 2 * N / 3, 0, N));
			print("x[%i] = ti2 + tr3;\n", index(o, k + 2 * N / 3, 1, N));
			break;
		case 4:
			print("const auto tr2 = zr1 + zr3;\n");
			print("const auto ti2 = zi1 + zi3;\n");
			print("const auto tr4 = zr1 - zr3;\n");
			print("const auto ti4 = zi1 - zi3;\n");
			print("x[%i] = zr0 + tr2;\n", index(o, k, 0, N));
			print("x[%i] = zi0 + ti2;\n", index(o, k, 1, N));
			print("x[%i] = zr2 + ti4;\n", index(o, k + N / 4, 0, N));
			print("x[%i] = zi2 - tr4;\n", index(o, k + N / 4, 1, N));
			print("x[%i] = zr0 - tr2;\n", index(o, k + N / 2, 0, N));
			print("x[%i] = zi0 - ti2;\n", index(o, k + N / 2, 1, N));
			print("x[%i] = zr2 - ti4;\n", index(o, k + 3 * N / 4, 0, N));
			print("x[%i] = zi2 + tr4;\n", index(o, k + 3 * N / 4, 1, N));
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
			print("x[%i] = zr0 + tr5;\n", index(o, k + 0 * N / 5, 0, N));
			print("x[%i] = zi0 + ti5;\n", index(o, k + 0 * N / 5, 1, N));
			print("x[%i] = tr8 + ti10;\n", index(o, k + 1 * N / 5, 0, N));
			print("x[%i] = ti8 - tr10;\n", index(o, k + 1 * N / 5, 1, N));
			print("x[%i] = tr9 + ti11;\n", index(o, k + 2 * N / 5, 0, N));
			print("x[%i] = ti9 - tr11;\n", index(o, k + 2 * N / 5, 1, N));
			print("x[%i] = tr9 - ti11;\n", index(o, k + 3 * N / 5, 0, N));
			print("x[%i] = ti9 + tr11;\n", index(o, k + 3 * N / 5, 1, N));
			print("x[%i] = tr8 - ti10;\n", index(o, k + 4 * N / 5, 0, N));
			print("x[%i] = ti8 + tr10;\n", index(o, k + 4 * N / 5, 1, N));
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
			print("x[%i] = tr7 + tr10;\n", index(o, k + 0 * N / 6, 0, N));
			print("x[%i] = tr8 - tr11;\n", index(o, k + 1 * N / 6, 0, N));
			print("x[%i] = tr9 + tr12;\n", index(o, k + 2 * N / 6, 0, N));
			print("x[%i] = tr7 - tr10;\n", index(o, k + 3 * N / 6, 0, N));
			print("x[%i] = tr8 + tr11;\n", index(o, k + 4 * N / 6, 0, N));
			print("x[%i] = tr9 - tr12;\n", index(o, k + 5 * N / 6, 0, N));
			print("x[%i] = ti7 + ti10;\n", index(o, k + 0 * N / 6, 1, N));
			print("x[%i] = ti8 - ti11;\n", index(o, k + 1 * N / 6, 1, N));
			print("x[%i] = ti9 + ti12;\n", index(o, k + 2 * N / 6, 1, N));
			print("x[%i] = ti7 - ti10;\n", index(o, k + 3 * N / 6, 1, N));
			print("x[%i] = ti8 + ti11;\n", index(o, k + 4 * N / 6, 1, N));
			print("x[%i] = ti9 - ti12;\n", index(o, k + 5 * N / 6, 1, N));
			break;

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
			print("x[%i] = ", index(o, k, 0, N));
			for (int i = 0; i < r; i++) {
				print_notab(" + zr%i", i);
			}
			print_notab(";\n");
			print("x[%i] = ", index(o, k, 1, N));
			for (int i = 0; i < r; i++) {
				print_notab(" + zi%i", i);
			}
			print_notab(";\n");
			for (int i = 1; i <= (r - 1) / 2; i++) {
				print("x[%i] = ap%i + am%i;\n", index(o, k + i * N / r, 0, N), i, i);
				print("x[%i] = bp%i - bm%i;\n", index(o, k + i * N / r, 1, N), i, i);
				print("x[%i] = ap%i - am%i;\n", index(o, k + (r - i) * N / r, 0, N), i, i);
				print("x[%i] = bp%i + bm%i;\n", index(o, k + (r - i) * N / r, 1, N), i, i);
			}
		}
		deindent();
		print("}\n");
	}
	deindent();
	print("}\n");
}

fft_type best_radix(int N, int o, bool first) {
	fft_type fftt;
	int best_cnt = 999999999;
	int best_radix = -1;
	for (int r = 2; r <= N; r++) {
		if (N % r == 0) {
			int this_cnt;
			if (r <= 6 || is_prime(r)) {
				this_cnt = fft_radix_opcnt(r, N);
				if (first) {
					this_cnt += N * MWEIGHT;
				}
				if (this_cnt < best_cnt) {
					best_cnt = this_cnt;
					fftt.N1 = r;
					fftt.type = RADIX;
					fftt.nops = best_cnt;
				}
			}
		}
	}
	if (N > 8 && is_prime(N)) {
		int raders_cnt = raders_fft_opcnt(N);
		if (raders_cnt < best_cnt) {
			best_cnt = raders_cnt;
			fftt.N1 = N;
			fftt.type = RADERS;
			fftt.nops = best_cnt;
		}
	}
	if (N > 6 ) {
		auto pfac = prime_fac(N);
		if (pfac.size() >= 2) {
			int N1 = 1, N2 = 1;
			int arrow = 1;
			for (auto i = pfac.begin(); i != pfac.end(); i++) {
				if (arrow > 0) {
					N1 *= pow(i->first, i->second);
				} else {
					N2 *= pow(i->first, i->second);
				}
				arrow = -arrow;
			}
			if (N1 > N2) {
				std::swap(N1, N2);
			}
			int gt_cnt = gt2_fft_opcnt(N1, N2);
			if (first) {
				gt_cnt += N * MWEIGHT;
			}
			if (gt_cnt < best_cnt) {
				fftt.type = GOOD;
				fftt.N1 = N1;
				fftt.N2 = N2;
				fftt.N3 = -1;
				fftt.nops = gt_cnt;
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
