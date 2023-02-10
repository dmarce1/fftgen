#include "fftgen.hpp"
#include <list>
#include <unordered_set>

std::vector<int> fft_bitreverse_indices_real(int N) {
	std::vector<int> I;
	std::vector<int> indices;
	I.resize(0);
	for (int n = 0; n < N; n++) {
		I.push_back(n);
	}
	indices = fft_bitr_real(N, 0, I, true);
	return indices;
}

std::vector<int> fft_bitr_real(int N, int o, std::vector<int> indices, bool first) {
	if (N == 1) {
		return indices;
	}
	auto fftt = best_radix_real(N, o, first);
	if (fftt.type == RADIX) {
		return fft_radix_bitr_real(fftt.N1, N, o, indices);
	} else if (fftt.type == GOOD) {
		return indices;
	} else {
		return indices;
	}
}

void gt2_fft_real(int N1, int N2, int o) {
	printf("// good-thomas real - %i = %i x %i\n", N1 * N2, N1, N2);
	print("{\n");
	indent();
	int N = N1 * N2;
	std::vector<std::vector<double>> y(N1, std::vector<double>(N2));
	std::vector<std::vector<std::complex<double>>>Y(N1, std::vector<std::complex<double>>(N2 / 2 + 1));
	std::vector<double> z(N1);
	std::vector<std::vector<std::complex<double>>>Z(N2 / 2 + 1, std::vector<std::complex<double>>(N1));
	std::vector<std::complex<double>> Z1(N1 / 2 + 1);
	print("std::array<double, %i> y;\n", N);
	print("std::array<double, %i> z0;\n", N1);
	print("std::array<std::array<double, %i>, %i> z1;\n", 2 * N1, N2 / 2);
	for (int n1 = 0; n1 < N1; n1++) {
		for (int n2 = 0; n2 < N2; n2++) {
//			y[n1][n2] = x[(N1 * n2 + N2 * n1) % N];
			print("y[%i] = x[%i];\n", n1 * N2 + n2, (N1 * n2 + N2 * n1) % N);
		}
	}
	for (int n1 = 0; n1 < N1; n1++) {
		print("{\n");
		indent();
		print("auto* x = y.data() + %i;", n1 * N2);
		fft_bitreverse_real(N2);
		fft_real(N2, 0);
		deindent();
		print("}\n");
	}
	for (int n1 = 0; n1 < N1; n1++) {
		print("z0[%i] = y[%i];\n", n1, n1 * N2);
		for (int n2 = 1; n2 < N2 / 2 + 1; n2++) {
			print("z1[%i][%i] = y[%i];\n", n2 - 1, 2 * n1, n1 * N2 + n2);
			if (n2 != 0 && !(N2 % 2 == 0 && n2 == N2 / 2)) {
				print("z1[%i][%i] = y[%i];\n", n2 - 1, 2 * n1 + 1, n1 * N2 + N2 - n2);
			}
		}
	}
	print("{\n");
	indent();
	print("auto* x = z0.data();\n");
	fft_bitreverse_real(N1);
	fft_real(N1, 0);
	deindent();
	print("}\n");
	for (int n2 = 1; n2 < N2 / 2 + 1; n2++) {
		print("{\n");
		indent();
		print("auto* x = z1[%i].data();\n", n2 - 1);
		fft_bitreverse(N1);
		fft(N1, 0);
		deindent();
		print("}\n");
	}
	for (int n = 0; n < N / 2 + 1; n++) {
		int n1 = n % N1;
		int n2 = n % N2;
		if (n2 == 0) {
			if (n1 >= N1 / 2 + 1) {
				print("x[%i] = z0[%i];\n", n, N1 - n1);
				if (n != 0 && !(N % 2 == 0 && n == N / 2)) {
					print("x[%i] = -z0[%i];\n", N - n, n1);
				}
			} else {
				print("x[%i] = z0[%i];\n", n, n1);
				if (n != 0 && !(N % 2 == 0 && n == N / 2)) {
					print("x[%i] = z0[%i];\n", N - n, N1 - n1);
				}
			}
		} else {
			if (n2 < N2 / 2 + 1) {
				print("x[%i] = z1[%i][%i];\n", n, n2 - 1, 2 * n1);
				if (n != 0 && !(N % 2 == 0 && n == N / 2)) {
					print("x[%i] = z1[%i][%i];\n", N - n, n2 - 1, 2 * n1 + 1);
				}
			} else {
				const int ni = n1 == 0 ? 0 : N1 - n1;
				print("x[%i] = z1[%i][%i];\n", n, N2 - n2 - 1, 2 * ni);
				if (n != 0 && !(N % 2 == 0 && n == N / 2)) {
					print("x[%i] = -z1[%i][%i];\n", N - n, N2 - n2 - 1, 2 * ni + 1);
				}

			}
		}
	}
	deindent();
	print("}\n");
}

std::vector<int> fft_radix_bitr_real(int r, int N, int o, std::vector<int> I) {
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
			auto K = fft_bitr_real(N / r, o + n * N / r, J);
			L.insert(L.end(), K.begin(), K.end());
		}
	};
	return L;
}

void fft_bitreverse_real(int N, std::vector<int> indices, int o) {
	if (indices.size() == 0) {
		indices = fft_bitreverse_indices_real(N);
	}
	std::vector<std::list<int>> strings;
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
				string.push_back(index_real(o, current, 0, N));
				while (first != next) {
					string.push_back(index_real(o, next, 0, N));
					current = next;
					touched[current] = true;
					next = indices[next];
				}
			}
			strings.push_back(string);
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

fft_type best_radix_real(int N, int o, bool first) {
	fft_type fftt;
	int best_cnt = 999999999;
	int best_radix = -1;
	for (int r = 2; r <= N; r++) {
		if (N % r == 0) {
			int this_cnt;
			if (r <= 6 || is_prime(r)) {
				this_cnt = fft_radix_real_opcnt(r, N, o);
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
	if (N >= 7 && is_prime(N)) {
		int raders_cnt = raders_fft_real_opcnt(N, 0);
		if (raders_cnt < best_cnt) {
			best_cnt = raders_cnt;
			fftt.N1 = N;
			fftt.type = RADERS;
			fftt.nops = best_cnt;
		}
	}
	if (N > 6) {
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
			//	int gt_cnt = gt2_fft_opcnt(N1, N2);
			if (first) {
				//			gt_cnt += N * MWEIGHT;
			}
			//		if (gt_cnt < best_cnt) {
			fftt.type = GOOD;
			fftt.N1 = N1;
			fftt.N2 = N2;
			fftt.N3 = -1;
			fftt.nops = 0;
			//	}
		}
	}
	return fftt;
}

void fft_real(int N, int o, bool first) {
	if (N == 1) {
		return;
	}
	auto fftt = best_radix_real(N, o, first);
	if (fftt.type == RADIX) {
		fft_radix_real(fftt.N1, N, o);
	} else if (fftt.type == RADERS) {
		raders_fft_real(fftt.N1, o);
	} else if (fftt.type == GOOD) {
		gt2_fft_real(fftt.N1, fftt.N2, o);
	}
}

int fft_real_opcnt(int N, int o) {
	fft_type fftt = best_radix_real(N, o);
	return fft_radix_real_opcnt(fftt.N1, N, o);
}

int mod(int n, int N) {
	return (n + N) % N;
}

void fft_complex_real(int N) {
	print("{\n");
	indent();
	print("std::array<double, %i> z;\n", 2 * (N / 2));
	fft_bitreverse(N / 2);
	fft(N / 2, 0, true);
	for (int k = 0; k <= N / 4; k++) {
		print("{\n");
		indent();
		print("const auto Xer = 0.5 * (x[%i] + x[%i]);\n", 2 * k, 2 * ((N / 2 - k) % (N / 2)));
		print("const auto Xei = 0.5 * (x[%i] - x[%i]);\n", 2 * k + 1, 2 * ((N / 2 - k) % (N / 2)) + 1);
		print("const auto Xoi = -0.5 * (x[%i] - x[%i]);\n", 2 * k, 2 * ((N / 2 - k) % (N / 2)));
		print("const auto Xor = 0.5 * (x[%i] + x[%i]);\n", 2 * k + 1, 2 * ((N / 2 - k) % (N / 2)) + 1);
		const auto W = twiddle(k, N);
		print("z[%i] = Xer + Xor * (%.17e) - Xoi * (%.17e);\n", k, W.real(), W.imag());
		if (k != 0) {
			print("z[%i] = Xer - Xor * (%.17e) + Xoi * (%.17e);\n", N / 2 - k, W.real(), W.imag());
			print("z[%i] = Xei + Xor * (%.17e) + Xoi * (%.17e);\n", N - k, W.imag(), W.real());
			print("z[%i] = -Xei + Xor * (%.17e) + Xoi * (%.17e);\n", N / 2 + k, W.imag(), W.real());
		} else {
			print("z[%i] = Xer - Xor;\n", N / 2);
		}
		deindent();
		print("}\n");
	}
	for (int k = 0; k < N; k++) {
		print("x[%i] = z[%i];\n", k, k);
	}
	deindent();
	print("}\n");
}

void fft_radix_real(int r, int N, int o, bool first) {
	if (N < r) {
		return;
	}
	print("{\n");
	indent();
	for (int n = 0; n < r; n++) {
		fft_real(N / r, o + n * N / r, first);
	}
	for (int k = 0; k < N / r; k++) {
		if (k > N / r / 2) {
			continue;
		}
		print("{\n");
		indent();
		std::vector<bool> rz(r, false);
		std::vector<bool> iz(r, false);
		if (k == 0 || (k == N / r / 2 && N / r % 2 == 0)) {
			for (int i = 0; i < r; i++) {
				const auto W = twiddle(i * k, N);
				const int L = N / r;
				const int M = L / 2;
				const int oo = o + i * N / r;
				const int iir = index_real(oo, k, 0, L);
				const int iii = index_real(oo, k, 1, L);
				if (i * k == 0) {
					print("const auto zr%i = x[%i];\n", i, iir);
					iz[i] = true;
				} else if (i * k == N / 8 && N % 8 == 0) {
					print("const auto zr%i = M_SQRT1_2 * (x[%i]);\n", i, iir);
					print("const auto zi%i = -zr%i;\n", i, i);
				} else if (i * k == N / 4 && N % 4 == 0) {
					print("const auto zi%i = -x[%i];\n", i, iir);
					rz[i] = true;
				} else if (i * k == 3 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = -M_SQRT1_2 * x[%i];\n", i, iir);
					print("const auto& zi%i = zr%i;\n", i, i);
				} else if (i * k == N / 2 && N % 2 == 0) {
					print("const auto zr%i = -x[%i];\n", i, iir);
					iz[i] = true;
				} else if (i * k == 5 * N / 8 && N % 8 == 0) {
					print("const auto zi%i = M_SQRT1_2 * x[%i];\n", i, iir);
					print("const auto zr%i = -zi%i;\n", i, i);
				} else if (i * k == 3 * N / 4 && N % 4 == 0) {
					print("const auto zi%i = x[%i];\n", i, iir);
					rz[i] = true;
				} else if (i * k == 7 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = M_SQRT1_2 * (x[%i]);\n", i, iir);
					print("const auto& zi%i = zr%i;\n", i, i);
				} else {
					print("const auto zr%i = (%.17e) * x[%i];\n", i, W.real(), iir);
					print("const auto zi%i = (%.17e) * x[%i];\n", i, W.imag(), iir);
				}
			}
		} else {
			const int M = N / 2;
			for (int i = 0; i < r; i++) {
				const auto W = twiddle(i * k, N);
				const int L = N / r;
				const int M = L / 2;
				const int oo = o + i * N / r;
				const int iir = index_real(oo, k, 0, L);
				const int iii = index_real(oo, k, 1, L);
				if (i * k == 0) {
					print("const auto zr%i = x[%i];\n", i, iir);
					print("const auto zi%i = x[%i];\n", i, iii);
				} else if (i * k == N / 8 && N % 8 == 0) {
					print("const auto zr%i = M_SQRT1_2 * (x[%i] + x[%i]);\n", i, iir, iii);
					print("const auto zi%i = -M_SQRT1_2 * (x[%i] - x[%i]);\n", i, iir, iii);
				} else if (i * k == N / 4 && N % 4 == 0) {
					print("const auto zr%i = x[%i];\n", i, iii);
					print("const auto zi%i = -x[%i];\n", i, iir);
				} else if (i * k == 3 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = -M_SQRT1_2 * (x[%i] - x[%i]);\n", i, iir, iii);
					print("const auto zi%i = -M_SQRT1_2 * (x[%i] + x[%i]);\n", i, iir, iii);
				} else if (i * k == N / 2 && N % 2 == 0) {
					print("const auto zr%i = -x[%i];\n", i, iir);
					print("const auto zi%i = -x[%i];\n", i, iii);
				} else if (i * k == 5 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = -M_SQRT1_2 * (x[%i] + x[%i]);\n", i, iir, iii);
					print("const auto zi%i = M_SQRT1_2 * (x[%i] - x[%i]);\n", i, iir, iii);
				} else if (i * k == 3 * N / 4 && N % 4 == 0) {
					print("const auto zr%i = -x[%i];\n", i, iii);
					print("const auto zi%i = x[%i];\n", i, iir);
				} else if (i * k == 7 * N / 8 && N % 8 == 0) {
					print("const auto zr%i = M_SQRT1_2 * (x[%i] - x[%i]);\n", i, iir, iii);
					print("const auto zi%i = M_SQRT1_2 * (x[%i] + x[%i]);\n", i, iir, iii);
				} else {
					print("const auto zr%i = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", i, W.real(), iir, -W.imag(), iii);
					print("const auto zi%i = std::fma((%.17e), x[%i], (%.17e) * x[%i]);\n", i, W.real(), iii, W.imag(), iir);
				}
			}
		}

		switch (r) {
		case 2: {
			print("// Radix - 2\n");
			std::unordered_set<int> used;
			int k0 = k;
			int k1 = N / 2 - k;
			if (used.find(k0) == used.end()) {
				used.insert(k0);
				if (rz[0] && rz[1]) {
				} else if (rz[0] && !rz[1]) {
					print("x[%i] = zr1;\n", index_real(o, k0, 0, N));
				} else if (!rz[0] && rz[1]) {
					print("x[%i] = zr0;\n", index_real(o, k0, 0, N));
				} else if (!rz[0] && !rz[1]) {
					print("x[%i] = zr0 + zr1;\n", index_real(o, k0, 0, N));
				}
			}
			if (used.find(k1) == used.end()) {
				used.insert(k1);
				if (rz[0] && rz[1]) {
				} else if (rz[0] && !rz[1]) {
					print("x[%i] = -zr1;\n", index_real(o, k1, 0, N));
				} else if (!rz[0] && rz[1]) {
					print("x[%i] = zr0;\n", index_real(o, k1, 0, N));
				} else if (!rz[0] && !rz[1]) {
					print("x[%i] = zr0 - zr1;\n", index_real(o, k1, 0, N));
				}
			}
			used.clear();
			if (k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					if (iz[0] && iz[1]) {
					} else if (iz[0] && !iz[1]) {
						print("x[%i] = zi1;\n", index_real(o, k0, 1, N));
					} else if (!iz[0] && iz[1]) {
						print("x[%i] = zi0;\n", index_real(o, k0, 1, N));
					} else if (!iz[0] && !iz[1]) {
						print("x[%i] = zi0 + zi1;\n", index_real(o, k0, 1, N));
					}
				}
			}
			if (k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					if (iz[0] && iz[1]) {
					} else if (!iz[0] && iz[1]) {
						print("x[%i] = -zi0;\n", index_real(o, k1, 1, N));
					} else if (iz[0] && iz[1]) {
						print("x[%i] = zi1;\n", index_real(o, k1, 1, N));
					} else if (!iz[0] && !iz[1]) {
						print("x[%i] = zi1 - zi0;\n", index_real(o, k1, 1, N));
					}
				}
			}
		}
			break;
		case 3:
			print("// Radix - 3\n");
			{
				bool tr1z = false;
				bool tr2z = false;
				bool tr3z = false;
				bool ti1z = false;
				bool ti2z = false;
				bool ti3z = false;
				if (rz[1] && rz[2]) {
					tr1z = true;
				} else if (!rz[1] && rz[2]) {
					print("const auto& tr1 = zr2;\n");
				} else if (rz[1] && !rz[2]) {
					print("const auto& tr1 = zr1;\n");
				} else {
					print("const auto tr1 = zr1 + zr2;\n");
				}
				if (tr1z && rz[0]) {
					tr2z = true;
				} else if (tr1z && !rz[0]) {
					print("const auto& tr2 = zr0;\n");
				} else if (!tr1z && rz[0]) {
					print("const auto tr2 = -0.5 * tr1;\n");
				} else {
					print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
				}
				if (rz[1] && rz[2]) {
					tr3z = true;
				} else if (!rz[1] && rz[2]) {
					print("const auto tr3 = (%24.17e) * zr2;\n", -sqrt(3) * 0.5);
				} else if (rz[1] && !rz[2]) {
					print("const auto tr3 = (%24.17e) * zr1;\n", sqrt(3) * 0.5);
				} else {
					print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
				}
				if (iz[1] && iz[2]) {
					ti1z = true;
				} else if (!iz[1] && iz[2]) {
					print("const auto& ti1 = zi2;\n");
				} else if (iz[1] && !iz[2]) {
					print("const auto& ti1 = zi1;\n");
				} else {
					print("const auto ti1 = zi1 + zi2;\n");
				}
				if (ti1z && iz[0]) {
					ti2z = true;
				} else if (ti1z && !iz[0]) {
					print("const auto& ti2 = zi0;\n");
				} else if (!ti1z && iz[0]) {
					print("const auto ti2 = -0.5 * ti1;\n");
				} else {
					print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
				}
				if (iz[1] && iz[2]) {
					ti3z = true;
				} else if (!iz[1] && iz[2]) {
					print("const auto ti3 = (%24.17e) * zi2;\n", -sqrt(3) * 0.5);
				} else if (iz[1] && !iz[2]) {
					print("const auto ti3 = (%24.17e) * zi1;\n", sqrt(3) * 0.5);
				} else {
					print("const auto ti3 = (%24.17e) * (zi1 - zi2);\n", sqrt(3) * 0.5);
				}
				const int k0 = k;
				const int k1 = k <= N / 2 ? k + N / 3 : N - (k + N / 3);
				const int k2 = N - (k + 2 * N / 3);
				std::unordered_set<int> used;
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					if (rz[0] && tr1z) {
					} else if (!rz[0] && tr1z) {
						print("x[%i] = zr0;\n", index_real(o, k0, 0, N));
					} else if (rz[0] && !tr1z) {
						print("x[%i] = tr1;\n", index_real(o, k0, 0, N));
					} else {
						print("x[%i] = zr0 + tr1;\n", index_real(o, k0, 0, N));
					}
				}
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					if (tr2z && ti3z) {
					} else if (!tr2z && ti3z) {
						print("x[%i] = tr2;\n", index_real(o, k1, 0, N));
					} else if (tr2z && !ti3z) {
						print("x[%i] = ti3;\n", index_real(o, k1, 0, N));
					} else {
						print("x[%i] = tr2 + ti3;\n", index_real(o, k1, 0, N));
					}
				}
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					if (tr2z && ti3z) {
					} else if (!tr2z && ti3z) {
						print("x[%i] = tr2;\n", index_real(o, k2, 0, N));
					} else if (tr2z && !ti3z) {
						print("x[%i] = -ti3;\n", index_real(o, k2, 0, N));
					} else {
						print("x[%i] = tr2 - ti3;\n", index_real(o, k2, 0, N));
					}
				}
				used.clear();
				if (k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
					if (used.find(k0) == used.end()) {
						used.insert(k0);
						if (iz[0] && ti1z) {
						} else if (!iz[0] && ti1z) {
							print("x[%i] = zi0;\n", index_real(o, k0, 1, N));
						} else if (iz[0] && !ti1z) {
							print("x[%i] = ti1;\n", index_real(o, k0, 1, N));
						} else {
							print("x[%i] = zi0 + ti1;\n", index_real(o, k0, 1, N));
						}
					}
				}
				if (k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					if (used.find(k1) == used.end()) {
						used.insert(k1);
						if (ti2z && tr3z) {
						} else if (!ti2z && tr3z) {
							print("x[%i] = ti2;\n", index_real(o, k1, 1, N));
						} else if (ti2z && !tr3z) {
							print("x[%i] = -tr3;\n", index_real(o, k1, 1, N));
						} else {
							print("x[%i] = ti2 - tr3;\n", index_real(o, k1, 1, N));
						}
					}
				}
				if (k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					if (used.find(k2) == used.end()) {
						used.insert(k2);
						if (ti2z && ti3z) {
						} else if (!ti2z && ti3z) {
							print("x[%i] = -ti2;\n", index_real(o, k2, 1, N));
						} else if (ti2z && !ti3z) {
							print("x[%i] = -tr3;\n", index_real(o, k2, 1, N));
						} else {
							print("x[%i] = -(ti2 + tr3);\n", index_real(o, k2, 1, N));
						}
					}
				}
			}
			break;
		case 4: {
			print("// Radix - 4\n");
			const int k0 = k;
			const int k1 = k + N / 4;
			const int k2 = N - (k + 2 * N / 4);
			const int k3 = N - (k + 3 * N / 4);
			bool tr1z = false;
			bool tr2z = false;
			bool tr3z = false;
			bool tr4z = false;
			bool ti1z = false;
			bool ti2z = false;
			bool ti3z = false;
			bool ti4z = false;

			if (rz[0] && rz[2]) {
				tr1z = true;
				tr3z = true;
			} else if (!rz[0] && rz[2]) {
				print("const auto& tr1 = zr0;\n");
				print("const auto& tr3 = zr0;\n");
			} else if (rz[0] && !rz[2]) {
				print("const auto& tr1 = zr2;\n");
				print("const auto tr3 = -zr2;\n");
			} else {
				print("const auto tr1 = zr0 + zr2;\n");
				print("const auto tr3 = zr0 - zr2;\n");
			}
			if (rz[1] && rz[3]) {
				tr2z = true;
				tr4z = true;
			} else if (!rz[1] && rz[3]) {
				print("const auto& tr2 = zr1;\n");
				print("const auto& tr4 = zr1;\n");
			} else if (rz[1] && !rz[3]) {
				print("const auto& tr2 = zr3;\n");
				print("const auto tr4 = -zr3;\n");
			} else {
				print("const auto tr2 = zr1 + zr3;\n");
				print("const auto tr4 = zr1 - zr3;\n");
			}
			if (iz[0] && iz[2]) {
				ti1z = true;
				ti3z = true;
			} else if (!iz[0] && iz[2]) {
				print("const auto& ti1 = zi0;\n");
				print("const auto& ti3 = zi2;\n");
			} else if (iz[0] && !iz[2]) {
				print("const auto& ti1 = zi2;\n");
				print("const auto ti3 = -zi2;\n");
			} else {
				print("const auto ti1 = zi0 + zi2;\n");
				print("const auto ti3 = zi0 - zi2;\n");
			}
			if (iz[1] && iz[3]) {
				ti2z = true;
				ti4z = true;
			} else if (!iz[1] && iz[3]) {
				print("const auto& ti2 = zi1;\n");
				print("const auto& ti4 = zi1;\n");
			} else if (iz[1] && !iz[3]) {
				print("const auto& ti2 = zi3;\n");
				print("const auto ti4 = -zi3;\n");
			} else {
				print("const auto ti2 = zi1 + zi3;\n");
				print("const auto ti4 = zi1 - zi3;\n");
			}
			std::unordered_set<int> used;
			if (tr1z && tr2z) {
			} else if (!tr1z && tr2z) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					print("x[%i] = tr1;\n", index_real(o, k0, 0, N));
				}
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					print("x[%i] = tr1;\n", index_real(o, k2, 0, N));
				}
			} else if (tr1z && !tr2z) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					print("x[%i] = tr2;\n", index_real(o, k0, 0, N));
				}
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					print("x[%i] = -tr2;\n", index_real(o, k2, 0, N));
				}
			} else {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					print("x[%i] = tr1 + tr2;\n", index_real(o, k0, 0, N));
				}
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					print("x[%i] = tr1 - tr2;\n", index_real(o, k2, 0, N));
				}
			}

			if (tr3z && ti4z) {
			} else if (!tr3z && ti4z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					print("x[%i] = tr3;\n", index_real(o, k1, 0, N));
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					print("x[%i] = tr3;\n", index_real(o, k3, 0, N));
				}
			} else if (tr3z && !ti4z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					print("x[%i] = ti4;\n", index_real(o, k1, 0, N));
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					print("x[%i] = -ti4;\n", index_real(o, k3, 0, N));
				}
			} else {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					print("x[%i] = tr3 + ti4;\n", index_real(o, k1, 0, N));
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					print("x[%i] = tr3 - ti4;\n", index_real(o, k3, 0, N));
				}
			}
			used.clear();
			if (ti1z && ti2z) {
			} else if (!ti1z && ti2z) {
				if (used.find(k0) == used.end() && k0 != N - k0 && k0 != 0) {
					print("x[%i] = ti1;\n", index_real(o, k0, 1, N));
					used.insert(k0);
				}
				if (used.find(k2) == used.end() && k2 != N - k2 && k2 != 0) {
					print("x[%i] = -ti1;\n", index_real(o, k2, 1, N));
					used.insert(k2);
				}
			} else if (ti1z && !ti2z) {
				if (used.find(k0) == used.end() && k0 != N - k0 && k0 != 0) {
					print("x[%i] = ti2;\n", index_real(o, k0, 1, N));
					used.insert(k0);
				}
				if (used.find(k2) == used.end() && k2 != N - k2 && k2 != 0) {
					print("x[%i] = ti2;\n", index_real(o, k2, 1, N));
					used.insert(k2);
				}
			} else {
				if (used.find(k0) == used.end() && k0 != N - k0 && k0 != 0) {
					print("x[%i] = ti1 + ti2;\n", index_real(o, k0, 1, N));
					used.insert(k0);
				}
				if (used.find(k2) == used.end() && k2 != N - k2 && k2 != 0) {
					print("x[%i] = -(ti1 - ti2);\n", index_real(o, k2, 1, N));
					used.insert(k2);
				}
			}
			if (ti3z && tr4z) {
			} else if (!ti3z && tr4z) {
				if (used.find(k1) == used.end() && k1 != N - k1 && k1 != 0) {
					print("x[%i] = ti3;\n", index_real(o, k1, 1, N));
					used.insert(k1);
				}
				if (used.find(k3) == used.end() && k3 != N - k3 && k3 != 0) {
					print("x[%i] = -ti3;\n", index_real(o, k3, 1, N));
					used.insert(k3);
				}
			} else if (ti3z && !tr4z) {
				if (used.find(k1) == used.end() && k1 != N - k1 && k1 != 0) {
					print("x[%i] = -tr4;\n", index_real(o, k1, 1, N));
					used.insert(k1);
				}
				if (used.find(k3) == used.end() && k3 != N - k3 && k3 != 0) {
					print("x[%i] = -tr4;\n", index_real(o, k3, 1, N));
					used.insert(k3);
				}
			} else {
				if (used.find(k1) == used.end() && k1 != N - k1 && k1 != 0) {
					print("x[%i] = ti3 - tr4;\n", index_real(o, k1, 1, N));
					used.insert(k1);
				}
				if (used.find(k3) == used.end() && k3 != N - k3 && k3 != 0) {
					print("x[%i] = -(ti3 + tr4);\n", index_real(o, k3, 1, N));
					used.insert(k3);
				}
			}
		}
			break;
		case 5: {
			// Radix -- 5
			print("// Radix - 5\n");
			bool tr1z = false;
			bool tr2z = false;
			bool tr3z = false;
			bool tr4z = false;
			bool tr5z = false;
			bool tr6z = false;
			bool tr7z = false;
			bool tr8z = false;
			bool tr9z = false;
			bool tr10z = false;
			bool tr11z = false;
			bool ti1z = false;
			bool ti2z = false;
			bool ti3z = false;
			bool ti4z = false;
			bool ti5z = false;
			bool ti6z = false;
			bool ti7z = false;
			bool ti8z = false;
			bool ti9z = false;
			bool ti10z = false;
			bool ti11z = false;
			if (rz[1] && rz[4]) {
				tr1z = tr3z = true;
			} else if (!rz[1] && rz[4]) {
				print("const auto& tr1 = zr1;\n");
				print("const auto& tr3 = zr1;\n");
			} else if (rz[1] && !rz[4]) {
				print("const auto& tr1 = zr4;\n");
				print("const auto tr3 = -zr4;\n");
			} else {
				print("const auto tr1 = zr1 + zr4;\n");
				print("const auto tr3 = zr1 - zr4;\n");
			}
			if (rz[2] && rz[3]) {
				tr2z = tr4z = true;
			} else if (!rz[2] && rz[3]) {
				print("const auto& tr2 = zr2;\n");
				print("const auto& tr4 = zr2;\n");
			} else if (rz[2] && !rz[3]) {
				print("const auto& tr2 = zr3;\n");
				print("const auto tr4 = -zr3;\n");
			} else {
				print("const auto tr2 = zr2 + zr3;\n");
				print("const auto tr4 = zr2 - zr3;\n");
			}
			if (tr1z && tr2z) {
				tr5z = true;
				tr6z = true;
			} else if (!tr1z && tr2z) {
				print("const auto& tr5 = tr1;\n");
				print("const auto tr6 = (%24.17e) * tr1;\n", sqrt(5) * 0.25);
			} else if (tr2z && !tr1z) {
				print("const auto& tr5 = tr2;\n");
				print("const auto tr6 = (%24.17e) * tr2;\n", sqrt(5) * 0.25);
			} else {
				print("const auto tr5 = tr1 + tr2;\n");
				print("const auto tr6 = (%24.17e) * (tr1 - tr2);\n", sqrt(5) * 0.25);
			}
			if (rz[0] && tr5z) {
				tr7z = true;
			} else if (!rz[0] && tr5z) {
				print("const auto& tr7 = zr0;\n");
			} else if (rz[0] && !tr5z) {
				print("const auto tr7 = -0.25 * tr5;\n");
			} else {
				print("const auto tr7 = std::fma(tr5, -0.25, zr0);\n");
			}
			if (tr6z && tr7z) {
				tr8z = true;
				tr9z = true;
			} else if (tr6z && !tr7z) {
				print("const auto& tr8 = tr7;\n");
				print("const auto& tr9 = tr7;\n");
			} else if (!tr6z && tr7z) {
				print("const auto& tr8 = tr6;\n");
				print("const auto tr9 = -tr6;\n");
			} else {
				print("const auto tr8 = tr7 + tr6;\n");
				print("const auto tr9 = tr7 - tr6;\n");
			}
			if (tr3z && tr4z) {
				tr10z = tr11z = true;
			} else if (!tr3z && tr4z) {
				print("const auto tr10 = (%24.17e) * tr3;\n", sin(2.0 * M_PI / 5.0));
				print("const auto tr11 = (%24.17e) * tr3;\n", sin(2.0 * M_PI / 10.0));
			} else if (tr3z && !tr4z) {
				print("const auto tr10 = (%24.17e) * tr4;\n", sin(2.0 * M_PI / 10.0));
				print("const auto tr11 = (%24.17e) * tr4;\n", -sin(2.0 * M_PI / 5.0));
			} else {
				print("const auto tr10 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
				print("const auto tr11 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
			}

			if (iz[1] && iz[4]) {
				ti1z = ti3z = true;
			} else if (!iz[1] && iz[4]) {
				print("const auto& ti1 = zi1;\n");
				print("const auto& ti3 = zi1;\n");
			} else if (iz[1] && !iz[4]) {
				print("const auto& ti1 = zi4;\n");
				print("const auto ti3 = -zi4;\n");
			} else {
				print("const auto ti1 = zi1 + zi4;\n");
				print("const auto ti3 = zi1 - zi4;\n");
			}
			if (iz[2] && iz[3]) {
				ti2z = ti4z = true;
			} else if (!iz[2] && iz[3]) {
				print("const auto& ti2 = zi2;\n");
				print("const auto& ti4 = zi2;\n");
			} else if (iz[2] && !iz[3]) {
				print("const auto& ti2 = zi3;\n");
				print("const auto ti4 = -zi3;\n");
			} else {
				print("const auto ti2 = zi2 + zi3;\n");
				print("const auto ti4 = zi2 - zi3;\n");
			}
			if (ti1z && ti2z) {
				ti5z = true;
				ti6z = true;
			} else if (!ti1z && ti2z) {
				print("const auto& ti5 = ti1;\n");
				print("const auto ti6 = (%24.17e) * ti1;\n", sqrt(5) * 0.25);
			} else if (ti2z && !ti1z) {
				print("const auto& ti5 = ti2;\n");
				print("const auto ti6 = (%24.17e) * -ti2;\n", -sqrt(5) * 0.25);
			} else {
				print("const auto ti5 = ti1 + ti2;\n");
				print("const auto ti6 = (%24.17e) * (ti1 - ti2);\n", sqrt(5) * 0.25);
			}
			if (iz[0] && ti5z) {
				ti7z = true;
			} else if (!iz[0] && ti5z) {
				print("const auto& ti7 = zi0;\n");
			} else if (iz[0] && !ti5z) {
				print("const auto ti7 = -0.25 * ti5;\n");
			} else {
				print("const auto ti7 = std::fma(ti5, -0.25, zi0);\n");
			}
			if (ti6z && ti7z) {
				ti8z = true;
				ti9z = true;
			} else if (ti6z && !ti7z) {
				print("const auto& ti8 = ti7;\n");
				print("const auto& ti9 = ti7;\n");
			} else if (!ti6z && ti7z) {
				print("const auto& ti8 = ti6;\n");
				print("const auto ti9 = -ti6;\n");
			} else {
				print("const auto ti8 = ti7 + ti6;\n");
				print("const auto ti9 = ti7 - ti6;\n");
			}
			if (ti3z && ti4z) {
				ti10z = ti11z = true;
			} else if (!ti3z && ti4z) {
				print("const auto ti10 = (%24.17e) * ti3;\n", sin(2.0 * M_PI / 5.0));
				print("const auto ti11 = (%24.17e) * ti3;\n", sin(2.0 * M_PI / 10.0));
			} else if (ti3z && !ti4z) {
				print("const auto ti10 = (%24.17e) * ti4;\n", sin(2.0 * M_PI / 10.0));
				print("const auto ti11 = (%24.17e) * ti4;\n", -sin(2.0 * M_PI / 5.0));
			} else {
				print("const auto ti10 = std::fma((%24.17e), ti3, (%24.17e) * ti4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
				print("const auto ti11 = std::fma((%24.17e), ti3, (%24.17e) * ti4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
			}
			int k0 = k;
			int k1 = k + N / 5;
			int k2 = k + 2 * N / 5;
			int k3 = N - (k + 3 * N / 5);
			int k4 = N - (k + 4 * N / 5);
			if (k > N / 2) {
				k2 = N - k2;
			}
			std::unordered_set<int> used;
			if (used.find(k0) == used.end()) {
				used.insert(k0);
				if (rz[0] && tr5z) {
				} else if (!rz[0] && tr5z) {
					print("x[%i] = zr0;\n", index_real(o, k0, 0, N));
				} else if (rz[0] && !tr5z) {
					print("x[%i] = tr5;\n", index_real(o, k0, 0, N));
				} else {
					print("x[%i] = zr0 + tr5;\n", index_real(o, k0, 0, N));
				}
			}
			if (tr8z && ti10z) {
			} else if (!tr8z && ti10z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					print("x[%i] = tr8;\n", index_real(o, k1, 0, N));
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
					print("x[%i] = tr8;\n", index_real(o, k4, 0, N));
				}
			} else if (tr8z && !ti10z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					print("x[%i] = ti10;\n", index_real(o, k1, 0, N));
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
					print("x[%i] = -ti10;\n", index_real(o, k4, 0, N));
				}
			} else {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					print("x[%i] = tr8 + ti10;\n", index_real(o, k1, 0, N));
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
					print("x[%i] = tr8 - ti10;\n", index_real(o, k4, 0, N));
				}
			}
			if (tr9z && ti11z) {
			} else if (!tr9z && ti11z) {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					print("x[%i] = tr9;\n", index_real(o, k2, 0, N));
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					print("x[%i] = tr9;\n", index_real(o, k3, 0, N));
				}
			} else if (tr9z && !ti11z) {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					print("x[%i] = ti11;\n", index_real(o, k2, 0, N));
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					print("x[%i] = -ti11;\n", index_real(o, k3, 0, N));
				}
			} else {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					print("x[%i] = tr9 + ti11;\n", index_real(o, k2, 0, N));
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					print("x[%i] = tr9 - ti11;\n", index_real(o, k3, 0, N));
				}
			}

			used.clear();
			if (used.find(k0) == used.end() && k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
				used.insert(k0);
				if (iz[0] && ti5z) {
				} else if (!iz[0] && ti5z) {
					print("x[%i] = zi0;\n", index_real(o, k0, 1, N));
				} else if (iz[0] && !ti5z) {
					print("x[%i] = ti5;\n", index_real(o, k0, 1, N));
				} else {
					print("x[%i] = zi0 + ti5;\n", index_real(o, k0, 1, N));
				}
			}
			if (ti8z && tr10z) {
			} else if (!ti8z && tr10z) {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
					print("x[%i] = ti8;\n", index_real(o, k1, 1, N));
					print("x[%i] = -ti8;\n", index_real(o, k4, 1, N));
				}
			} else if (ti8z && !tr10z) {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
					print("x[%i] = -tr10;\n", index_real(o, k1, 1, N));
					print("x[%i] = -tr10;\n", index_real(o, k4, 1, N));
				}
			} else {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
					print("x[%i] = ti8 - tr10;\n", index_real(o, k1, 1, N));
					print("x[%i] = -(ti8 + tr10);\n", index_real(o, k4, 1, N));
				}
			}
			if (ti9z && tr11z) {
			} else if (!ti9z && tr11z) {
				if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					used.insert(k2);
					print("x[%i] = ti9;\n", index_real(o, k2, 1, N));
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					print("x[%i] = -ti9;\n", index_real(o, k3, 1, N));
				}
			} else if (ti9z && !tr11z) {
				if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					used.insert(k2);
					print("x[%i] = -tr11;\n", index_real(o, k2, 1, N));
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					print("x[%i] = -tr11;\n", index_real(o, k3, 1, N));
				}
			} else {
				if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					used.insert(k2);
					print("x[%i] = ti9 - tr11;\n", index_real(o, k2, 1, N));
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					print("x[%i] = -(ti9 + tr11);\n", index_real(o, k3, 1, N));
				}
			}
		}
			break;
		case 6: {
			print("// Radix - 6\n");
			bool tr1z = false;
			bool tr2z = false;
			bool tr3z = false;
			bool tr4z = false;
			bool tr5z = false;
			bool tr6z = false;
			bool tr7z = false;
			bool tr8z = false;
			bool tr9z = false;
			bool tr10z = false;
			bool tr11z = false;
			bool tr12z = false;
			bool ti1z = false;
			bool ti2z = false;
			bool ti3z = false;
			bool ti4z = false;
			bool ti5z = false;
			bool ti6z = false;
			bool ti7z = false;
			bool ti8z = false;
			bool ti9z = false;
			bool ti10z = false;
			bool ti11z = false;
			bool ti12z = false;

			if (rz[2] && rz[4]) {
				tr1z = true;
				tr3z = true;
			} else if (!rz[2] && rz[4]) {
				print("const auto tr1 = zr2;\n");
				print("const auto tr3 = (%24.17e) * zr2;\n", sin(M_PI / 3.0));
			} else if (rz[2] && !rz[4]) {
				print("const auto tr1 = zr4;\n");
				print("const auto tr3 = (%24.17e) * zr4;\n", -sin(M_PI / 3.0));
			} else {
				print("const auto tr1 = zr2 + zr4;\n");
				print("const auto tr3 = (%24.17e) * (zr2 - zr4);\n", sin(M_PI / 3.0));
			}
			if (tr1z && rz[0]) {
				tr2z = true;
			} else if (!tr1z && rz[0]) {
				print("const auto tr2 = -0.5 * tr1;\n");
			} else if (tr1z && !rz[0]) {
				print("const auto& tr2 = zr0;\n");
			} else {
				print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
			}
			if (rz[1] && rz[5]) {
				tr4z = true;
				tr6z = true;
			} else if (!rz[1] && rz[5]) {
				print("const auto& tr4 = zr1;\n");
				print("const auto tr6 = (%24.17e) * zr1;\n", -sin(M_PI / 3.0));
			} else if (rz[1] && !rz[5]) {
				print("const auto& tr4 = zr5;\n");
				print("const auto tr6 = (%24.17e) * zr5;\n", sin(M_PI / 3.0));
			} else {
				print("const auto tr4 = zr5 + zr1;\n");
				print("const auto tr6 = (%24.17e) * (zr5 - zr1);\n", sin(M_PI / 3.0));
			}
			if (tr4z && rz[3]) {
				tr5z = true;
			} else if (!tr4z && rz[3]) {
				print("const auto tr5 = -0.5 * tr4;\n");
			} else if (tr4z && !rz[3]) {
				print("const auto& tr5 = zr3;\n");
			} else {
				print("const auto tr5 = std::fma(tr4, -0.5, zr3);\n");
			}

			if (iz[2] && iz[4]) {
				ti1z = true;
				ti3z = true;
			} else if (!iz[2] && iz[4]) {
				print("const auto ti1 = zi2;\n");
				print("const auto ti3 = (%24.17e) * zi2;\n", sin(M_PI / 3.0));
			} else if (iz[2] && !iz[4]) {
				print("const auto ti1 = zi4;\n");
				print("const auto ti3 = (%24.17e) * zi4;\n", -sin(M_PI / 3.0));
			} else {
				print("const auto ti1 = zi2 + zi4;\n");
				print("const auto ti3 = (%24.17e) * (zi2 - zi4);\n", sin(M_PI / 3.0));
			}
			if (ti1z && iz[0]) {
				ti2z = true;
			} else if (!ti1z && iz[0]) {
				print("const auto ti2 = -0.5 * ti1;\n");
			} else if (ti1z && !iz[0]) {
				print("const auto& ti2 = zi0;\n");
			} else {
				print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
			}
			if (iz[1] && iz[5]) {
				ti4z = true;
				ti6z = true;
			} else if (!iz[1] && iz[5]) {
				print("const auto& ti4 = zi1;\n");
				print("const auto ti6 = (%24.17e) * zi1;\n", -sin(M_PI / 3.0));
			} else if (iz[1] && !iz[5]) {
				print("const auto& ti4 = zi5;\n");
				print("const auto ti6 = (%24.17e) * zi5;\n", sin(M_PI / 3.0));
			} else {
				print("const auto ti4 = zi5 + zi1;\n");
				print("const auto ti6 = (%24.17e) * (zi5 - zi1);\n", sin(M_PI / 3.0));
			}
			if (ti4z && iz[3]) {
				ti5z = true;
			} else if (!ti4z && iz[3]) {
				print("const auto ti5 = -0.5 * ti4;\n");
			} else if (tr4z && !rz[3]) {
				print("const auto& ti5 = zi3;\n");
			} else {
				print("const auto ti5 = std::fma(ti4, -0.5, zi3);\n");
			}

			if (rz[0] && tr1z) {
				tr7z = true;
			} else if (!rz[0] && tr1z) {
				print("const auto& tr7 = zr0;\n");
			} else if (rz[0] && !tr1z) {
				print("const auto& tr7 = tr1;\n");
			} else {
				print("const auto tr7 = zr0 + tr1;\n");
			}

			if (iz[0] && ti1z) {
				ti7z = true;
			} else if (!iz[0] && ti1z) {
				print("const auto& ti7 = zi0;\n");
			} else if (iz[0] && !ti1z) {
				print("const auto& ti7 = ti1;\n");
			} else {
				print("const auto ti7 = zi0 + ti1;\n");
			}

			if (tr2z && ti3z) {
				tr8z = true;
				tr9z = true;
			} else if (!tr2z && ti3z) {
				print("const auto& tr8 = tr2;\n");
				print("const auto& tr9 = tr2;\n");
			} else if (tr2z && !ti3z) {
				print("const auto& tr8 = ti3;\n");
				print("const auto tr9 = -ti3;\n");
			} else {
				print("const auto tr8 = tr2 + ti3;\n");
				print("const auto tr9 = tr2 - ti3;\n");
			}

			if (ti2z && tr3z) {
				ti2z = true;
				tr3z = true;
			} else if (!ti2z && tr3z) {
				print("const auto& ti8 = ti2;\n");
				print("const auto& ti9 = ti2;\n");
			} else if (ti2z && !tr3z) {
				print("const auto ti8 = -tr3;\n");
				print("const auto& ti9 = tr3;\n");
			} else {
				print("const auto ti8 = ti2 - tr3;\n");
				print("const auto ti9 = ti2 + tr3;\n");
			}

			if (rz[3] && tr4z) {
				tr10z = true;
			} else if (rz[3] && !tr4z) {
				print("const auto& tr10 = tr4;\n");
			} else if (!rz[3] && tr4z) {
				print("const auto& tr10 = zr3;\n");
			} else {
				print("const auto tr10 = zr3 + tr4;\n");
			}

			if (iz[3] && ti4z) {
				ti4z = true;
			} else if (iz[3] && !ti4z) {
				print("const auto& ti10 = ti4;\n");
			} else if (!iz[3] && ti4z) {
				print("const auto& ti10 = zi3;\n");
			} else {
				print("const auto ti10 = zi3 + ti4;\n");
			}

			if (tr5z && ti6z) {
				tr5z = true;
				ti6z = true;
			} else if (!tr5z && ti6z) {
				print("const auto& tr11 = tr5;\n");
				print("const auto& tr12 = tr5;\n");
			} else if (tr5z && !ti6z) {
				print("const auto& tr11 = ti6;\n");
				print("const auto tr12 = -ti6;\n");
			} else {
				print("const auto tr11 = tr5 + ti6;\n");
				print("const auto tr12 = tr5 - ti6;\n");
			}

			if (ti5z && tr6z) {
				ti5z = true;
				tr6z = true;
			} else if (!ti5z && tr6z) {
				print("const auto& ti11 = ti5;\n");
				print("const auto& ti12 = ti5;\n");
			} else if (ti5z && !tr6z) {
				print("const auto ti11 = -tr6;\n");
				print("const auto& ti12 = tr6;\n");
			} else {
				print("const auto ti11 = ti5 - tr6;\n");
				print("const auto ti12 = ti5 + tr6;\n");
			}

			int k0 = k + 0 * N / 6;
			int k1 = k + 1 * N / 6;
			int k2 = k + 2 * N / 6;
			int k3 = N - (k + 3 * N / 6);
			int k4 = N - (k + 4 * N / 6);
			int k5 = N - (k + 5 * N / 6);
			std::unordered_set<int> used;
			if (tr7z && tr10z) {
			} else if (!tr7z && tr10z) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					print("x[%i] = tr7;\n", index_real(o, k0, 0, N));
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					print("x[%i] = tr7;\n", index_real(o, k3, 0, N));
				}
			} else if (tr7z && !tr10z) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					print("x[%i] = tr10;\n", index_real(o, k0, 0, N));
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					print("x[%i] = -tr10;\n", index_real(o, k3, 0, N));
				}
			} else {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					print("x[%i] = tr7 + tr10;\n", index_real(o, k0, 0, N));
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					print("x[%i] = tr7 - tr10;\n", index_real(o, k3, 0, N));
				}
			}
			if (tr8z && tr11z) {

			} else if (!tr8z && tr11z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					print("x[%i] = tr8;\n", index_real(o, k1, 0, N));
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
					print("x[%i] = tr8;\n", index_real(o, k4, 0, N));
				}
			} else if (tr8z && !tr11z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					print("x[%i] = -tr11;\n", index_real(o, k1, 0, N));
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
					print("x[%i] = tr11;\n", index_real(o, k4, 0, N));
				}
			} else {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					print("x[%i] = tr8 - tr11;\n", index_real(o, k1, 0, N));
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
					print("x[%i] = tr8 + tr11;\n", index_real(o, k4, 0, N));
				}
			}
			if (tr9z && tr12z) {
			} else if (tr9z && !tr12z) {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					print("x[%i] = tr12;\n", index_real(o, k2, 0, N));
				}
				if (used.find(k5) == used.end()) {
					used.insert(k5);
					print("x[%i] = -tr12;\n", index_real(o, k5, 0, N));
				}
			} else if (!tr9z && tr12z) {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					print("x[%i] = tr9;\n", index_real(o, k2, 0, N));
				}
				if (used.find(k5) == used.end()) {
					used.insert(k5);
					print("x[%i] = tr9;\n", index_real(o, k5, 0, N));
				}
			} else {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					print("x[%i] = tr9 + tr12;\n", index_real(o, k2, 0, N));
				}
				if (used.find(k5) == used.end()) {
					used.insert(k5);
					print("x[%i] = tr9 - tr12;\n", index_real(o, k5, 0, N));
				}
			}
			used.clear();

			if (ti7z && ti10z) {
			} else if (!ti7z && ti10z) {
				if (used.find(k0) == used.end() && k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
					used.insert(k0);
					print("x[%i] = ti7;\n", index_real(o, k0, 1, N));
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					print("x[%i] = -ti7;\n", index_real(o, k3, 1, N));
				}
			} else if (ti7z && !ti10z) {
				if (used.find(k0) == used.end() && k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
					used.insert(k0);
					print("x[%i] = ti10;\n", index_real(o, k0, 1, N));
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					print("x[%i] = ti10;\n", index_real(o, k3, 1, N));
				}
			} else {
				if (used.find(k0) == used.end() && k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
					used.insert(k0);
					print("x[%i] = ti7 + ti10;\n", index_real(o, k0, 1, N));
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					print("x[%i] = ti10 - ti7;\n", index_real(o, k3, 1, N));
				}
			}

			if (ti8z && ti11z) {

			} else if (!ti8z && ti11z) {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
					print("x[%i] = ti8;\n", index_real(o, k1, 1, N));
				}
				if (used.find(k4) == used.end() && k4 != 0 && !(k4 == N / 2 && N % 2 == 0)) {
					used.insert(k4);
					print("x[%i] = -ti8;\n", index_real(o, k4, 1, N));
				}
			} else if (ti8z && !ti11z) {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
					print("x[%i] = -ti11;\n", index_real(o, k1, 1, N));
				}
				if (used.find(k4) == used.end() && k4 != 0 && !(k4 == N / 2 && N % 2 == 0)) {
					used.insert(k4);
					print("x[%i] = -ti11;\n", index_real(o, k4, 1, N));
				}
			} else {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
					print("x[%i] = ti8 - ti11;\n", index_real(o, k1, 1, N));
				}
				if (used.find(k4) == used.end() && k4 != 0 && !(k4 == N / 2 && N % 2 == 0)) {
					used.insert(k4);
					print("x[%i] = -ti8 - ti11;\n", index_real(o, k4, 1, N));
				}
			}
			if (ti9z && ti12z) {
			} else if (!ti9z && ti12z) {
				if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					used.insert(k2);
					print("x[%i] = ti9;\n", index_real(o, k2, 1, N));
				}
				if (used.find(k5) == used.end() && k5 != 0 && !(k5 == N / 2 && N % 2 == 0)) {
					used.insert(k5);
					print("x[%i] = -ti9;\n", index_real(o, k5, 1, N));
				}
			} else if (ti9z && !ti12z) {
				if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					used.insert(k2);
					print("x[%i] = ti12;\n", index_real(o, k2, 1, N));
				}
				if (used.find(k5) == used.end() && k5 != 0 && !(k5 == N / 2 && N % 2 == 0)) {
					used.insert(k5);
					print("x[%i] = ti12;\n", index_real(o, k5, 1, N));
				}
			} else {
				if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					used.insert(k2);
					print("x[%i] = ti9 + ti12;\n", index_real(o, k2, 1, N));
				}
				if (used.find(k5) == used.end() && k5 != 0 && !(k5 == N / 2 && N % 2 == 0)) {
					used.insert(k5);
					print("x[%i] = ti12 - ti9;\n", index_real(o, k5, 1, N));
				}
			}

		}
			break;
		default: {
			std::vector<bool> txpz((r - 1) / 2 + 1, false);
			std::vector<bool> txmz((r - 1) / 2 + 1, false);
			std::vector<bool> typz((r - 1) / 2 + 1, false);
			std::vector<bool> tymz((r - 1) / 2 + 1, false);
			print("// radix - %i\n", r);
			for (int j = 1; j <= (r - 1) / 2; j++) {
				if (rz[j] && rz[r - j]) {
					txpz[j] = true;
				} else if (!rz[j] && rz[r - j]) {
					print("const auto txp%i = zr%i;\n", j, j);
				} else if (rz[j] && !rz[r - j]) {
					print("const auto txp%i = zr%i;\n", j, r - j);
				} else {
					print("const auto txp%i = zr%i + zr%i;\n", j, j, r - j);
				}
			}
			for (int j = 1; j <= (r - 1) / 2; j++) {
				if (rz[j] && rz[r - j]) {
					txmz[j] = true;
				} else if (!rz[j] && rz[r - j]) {
					print("const auto txm%i = zr%i;\n", j, j);
				} else if (rz[j] && !rz[r - j]) {
					print("const auto txm%i = -zr%i;\n", j, r - j);
				} else {
					print("const auto txm%i = zr%i - zr%i;\n", j, j, r - j);
				}
			}
			for (int j = 1; j <= (r - 1) / 2; j++) {
				if (iz[j] && iz[r - j]) {
					typz[j] = true;
				} else if (!iz[j] && iz[r - j]) {
					print("const auto typ%i = zi%i;\n", j, j);
				} else if (iz[j] && !iz[r - j]) {
					print("const auto typ%i = zi%i;\n", j, r - j);
				} else {
					print("const auto typ%i = zi%i + zi%i;\n", j, j, r - j);
				}
			}
			for (int j = 1; j <= (r - 1) / 2; j++) {
				if (iz[j] && iz[r - j]) {
					tymz[j] = true;
				} else if (!iz[j] && iz[r - j]) {
					print("const auto tym%i = zi%i;\n", j, j);
				} else if (iz[j] && !iz[r - j]) {
					print("const auto tym%i = -zi%i;\n", j, r - j);
				} else {
					print("const auto tym%i = zi%i - zi%i;\n", j, j, r - j);
				}
			}
			std::vector<bool> apz((r - 1) / 2 + 1, true);
			std::vector<bool> amz((r - 1) / 2 + 1, true);
			std::vector<bool> bpz((r - 1) / 2 + 1, true);
			std::vector<bool> bmz((r - 1) / 2 + 1, true);
			for (int i = 1; i <= (r - 1) / 2; i++) {
				bool apfirst = true;
				bool amfirst = true;
				bool bpfirst = true;
				bool bmfirst = true;

				for (int j = 1; j <= (r - 1) / 2; j++) {
					if (!txpz[j]) {
						apz[i] = false;
						if (apfirst) {
							apfirst = false;
							print("auto ap%i = txp%i * (%24.17e);\n", i, j, cos(2.0 * M_PI * j * i / r));
						} else {
							print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / r), i);
						}
					}
					if (!typz[j]) {
						bpz[i] = false;
						if (bpfirst) {
							bpfirst = false;
							print("auto bp%i = typ%i * (%24.17e);\n", i, j, cos(2.0 * M_PI * j * i / r));
						} else {
							print("bp%i = std::fma(typ%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / r), i);
						}
					}
					if (!tymz[j]) {
						amz[i] = false;
						if (amfirst) {
							print("auto am%i = tym%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
							amfirst = false;
						} else {
							print("am%i = std::fma(tym%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
						}
					}
					if (!txmz[j]) {
						bmz[i] = false;
						if (bmfirst) {
							print("auto bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
							bmfirst = false;
						} else {
							print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
						}
					}
				}
				if (!rz[0]) {
					apz[i] = false;
					if (apfirst) {
						apfirst = false;
						print("auto ap%i = zr0;\n", i);
					} else {
						print("ap%i += zr0;\n", i);
					}
				}
				if (!iz[0]) {
					bpz[i] = false;
					if (bpfirst) {
						bpfirst = false;
						print("auto bp%i = zi0;\n", i);
					} else {
						print("bp%i += zi0;\n", i);
					}
				}
			}
			std::unordered_set<int> used;
			if (used.find(k) == used.end()) {
				used.insert(k);
				print("x[%i] = ", index_real(o, k, 0, N));
				for (int i = 0; i < r; i++) {
					if (!rz[i]) {
						print_notab(" + zr%i", i);
					}
				}
				print_notab(";\n");
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				const int ka = k + i * N / r;
				const int rk0 = index_real(o, ka, 0, N);
				if (used.find(ka) == used.end()) {
					used.insert(ka);
					if (apz[i] && amz[i]) {
					} else if (apz[i] && !amz[i]) {
						print("x[%i] = am%i;i\n", rk0, i);
					} else if (!apz[i] && amz[i]) {
						print("x[%i] = ap%i; \n", rk0, i);
					} else {
						print("x[%i] = ap%i + am%i;\n", rk0, i, i);
					}
				}
			}
			if (!(k == N / r / 2 && N % 2 == 0)) {
				for (int i = 1; i <= (r - 1) / 2; i++) {
					const int kb = N - (k + (r - i) * N / r);
					const int rk1 = index_real(o, kb, 0, N);
					if (used.find(kb) == used.end()) {
						used.insert(kb);
						if (apz[i] && amz[i]) {
						} else if (apz[i] && !amz[i]) {
							print("x[%i] = -am%i; \n", rk1, i);
						} else if (!apz[i] && amz[i]) {
							print("x[%i] = ap%i; \n", rk1, i);
						} else {
							print("x[%i] = ap%i - am%i;\n", rk1, i, i);
						}
					}
				}
			}
			used.clear();
			if (used.find(k) == used.end()) {
				if (k != 0 && !(k == N / 2 && N % 2 == 0)) {
					used.insert(k);
					print("x[%i] = ", index_real(o, k, 1, N));
					for (int i = 0; i < r; i++) {
						if (!iz[i]) {
							print_notab(" + zi%i", i);
						}
					}
					print_notab(";\n");
				}
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				const int ka = k + i * N / r;
				const int ik0 = index_real(o, ka, 1, N);
				if (used.find(ka) == used.end() && ka != 0 && !(ka == N / 2 && N % 2 == 0)) {
					used.insert(ka);
					if (bpz[i] && bmz[i]) {
					} else if (bpz[i] && !bmz[i]) {
						print("x[%i] = -bm%i; \n", ik0, i);
					} else if (!bpz[i] && bmz[i]) {
						print("x[%i] = bp%i;  \n", ik0, i);
					} else {
						print("x[%i] = bp%i - bm%i; \n", ik0, i, i);
					}
				}
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				const int kb = N - (k + (r - i) * N / r);
				const int ik1 = index_real(o, kb, 1, N);
				if (used.find(kb) == used.end() && kb != 0 && !(kb == N / 2 && N % 2 == 0)) {
					used.insert(kb);
					if (bpz[i] && bmz[i]) {
					} else if (bpz[i] && !bmz[i]) {
						print("x[%i] = -bm%i; \n", ik1, i);
					} else if (!bpz[i] && bmz[i]) {
						print("x[%i] = -bp%i; \n", ik1, i);
					} else {
						print("x[%i] = -(bp%i + bm%i); \n", ik1, i, i);
					}
				}
			}
		}
		}
		deindent();
		print("}\n");
	}

	deindent();
	print("}\n");

}

int fft_radix_real_opcnt(int r, int N, int o) {
	int cnt = 0;
	if (N < r) {
		return cnt;
	}
	if (N > r) {
		for (int n = 0; n < r; n++) {
			cnt += fft_real_opcnt(N / r, o + n * N / r);
		}
	}
	for (int k = 0; k < N / r; k++) {
		if (k > N / r / 2) {
			continue;
		}
		std::vector<bool> rz(r, false);
		std::vector<bool> iz(r, false);
		if (k == 0 || (k == N / r / 2 && N / r % 2 == 0)) {
			for (int i = 0; i < r; i++) {
				const auto W = twiddle(i * k, N);
				const int L = N / r;
				const int M = L / 2;
				const int oo = o + i * N / r;
				const int iir = index_real(oo, k, 0, L);
				const int iii = index_real(oo, k, 1, L);
				if (i * k == 0) {
					iz[i] = true;
				} else if (i * k == N / 8 && N % 8 == 0) {
					cnt += 2;
				} else if (i * k == N / 4 && N % 4 == 0) {
					cnt++;
					rz[i] = true;
				} else if (i * k == 3 * N / 8 && N % 8 == 0) {
					cnt++;
				} else if (i * k == N / 2 && N % 2 == 0) {
					cnt++;
					iz[i] = true;
				} else if (i * k == 5 * N / 8 && N % 8 == 0) {
					cnt += 2;
				} else if (i * k == 3 * N / 4 && N % 4 == 0) {
					rz[i] = true;
				} else if (i * k == 7 * N / 8 && N % 8 == 0) {
					cnt++;
				} else {
					cnt += 2;
				}
			}
		} else {
			const int M = N / 2;
			for (int i = 0; i < r; i++) {
				const auto W = twiddle(i * k, N);
				const int L = N / r;
				const int M = L / 2;
				const int oo = o + i * N / r;
				const int iir = index_real(oo, k, 0, L);
				const int iii = index_real(oo, k, 1, L);
				if (i * k == 0) {
				} else if (i * k == N / 4 && N % 4 == 0) {
					cnt++;
				} else if (i * k == N / 2 && N % 2 == 0) {
					cnt += 2;
				} else if (i * k == 3 * N / 4 && N % 4 == 0) {
					cnt++;
				} else {
					cnt += 4;
				}
			}
		}

		switch (r) {
		case 2: {
			std::unordered_set<int> used;
			int k0 = k;
			int k1 = N / 2 - k;
			if (used.find(k0) == used.end()) {
				used.insert(k0);
				if (rz[0] && rz[1]) {
				} else if (rz[0] && !rz[1]) {
				} else if (!rz[0] && rz[1]) {
				} else if (!rz[0] && !rz[1]) {
					cnt++;
				}
			}
			if (used.find(k1) == used.end()) {
				used.insert(k1);
				if (rz[0] && rz[1]) {
				} else if (rz[0] && !rz[1]) {
					cnt++;
				} else if (!rz[0] && rz[1]) {
				} else if (!rz[0] && !rz[1]) {
					cnt++;
				}
			}
			used.clear();
			if (k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					if (iz[0] && iz[1]) {
					} else if (!iz[0] && iz[1]) {
					} else if (!iz[0] && !iz[1]) {
						cnt++;
					}
				}
			}
			if (k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					if (iz[0] && iz[1]) {
					} else if (!iz[0] && iz[1]) {
						cnt++;
					} else if (iz[0] && iz[1]) {
					} else if (!iz[0] && !iz[1]) {
						cnt++;
					}
				}
			}
		}
			break;
		case 3: {
			bool tr1z = false;
			bool tr2z = false;
			bool tr3z = false;
			bool ti1z = false;
			bool ti2z = false;
			bool ti3z = false;
			if (rz[1] && rz[2]) {
				tr1z = true;
			} else if (!rz[1] && rz[2]) {
			} else if (rz[1] && !rz[2]) {
			} else {
				cnt++;
			}
			if (tr1z && rz[0]) {
				tr2z = true;
			} else if (tr1z && !rz[0]) {
			} else if (!tr1z && rz[0]) {
				cnt++;
			} else {
				cnt++;
			}
			if (rz[1] && rz[2]) {
				tr3z = true;
			} else if (!rz[1] && rz[2]) {
				cnt++;
			} else if (rz[1] && !rz[2]) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (iz[1] && iz[2]) {
				ti1z = true;
			} else if (!iz[1] && iz[2]) {
			} else if (iz[1] && !iz[2]) {
			} else {
				cnt++;
			}
			if (ti1z && iz[0]) {
				ti2z = true;
			} else if (ti1z && !iz[0]) {
			} else if (!ti1z && iz[0]) {
				cnt++;
			} else {
				cnt++;
			}
			if (iz[1] && iz[2]) {
				ti3z = true;
			} else if (!iz[1] && iz[2]) {
				cnt++;
			} else if (iz[1] && !iz[2]) {
				cnt++;
			} else {
				cnt += 2;
			}
			const int k0 = k;
			const int k1 = k <= N / 2 ? k + N / 3 : N - (k + N / 3);
			const int k2 = N - (k + 2 * N / 3);
			std::unordered_set<int> used;
			if (used.find(k0) == used.end()) {
				used.insert(k0);
				if (rz[0] && tr1z) {
				} else if (!rz[0] && tr1z) {
				} else if (rz[0] && !tr1z) {
				} else {
					cnt++;
				}
			}
			if (used.find(k1) == used.end()) {
				used.insert(k1);
				if (tr2z && ti3z) {
				} else if (!tr2z && ti3z) {
				} else if (tr2z && !ti3z) {
				} else {
					cnt++;
				}
			}
			if (used.find(k2) == used.end()) {
				used.insert(k2);
				if (tr2z && ti3z) {
				} else if (!tr2z && ti3z) {
				} else if (tr2z && !ti3z) {
					cnt++;
				} else {
					cnt++;
				}
			}
			used.clear();
			if (k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					if (iz[0] && ti1z) {
					} else if (!iz[0] && ti1z) {
					} else if (iz[0] && !ti1z) {
					} else {
						cnt++;
					}
				}
			}
			if (k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					if (ti2z && tr3z) {
					} else if (!ti2z && tr3z) {
					} else if (ti2z && !tr3z) {
						cnt++;
					} else {
						cnt++;
					}
				}
			}
			if (k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					if (ti2z && ti3z) {
					} else if (!ti2z && ti3z) {
						cnt++;
					} else if (ti2z && !ti3z) {
						cnt++;
					} else {
						cnt += 2;
					}
				}
			}
		}
			break;
		case 4: {
			const int k0 = k;
			const int k1 = k + N / 4;
			const int k2 = N - (k + 2 * N / 4);
			const int k3 = N - (k + 3 * N / 4);
			bool tr1z = false;
			bool tr2z = false;
			bool tr3z = false;
			bool tr4z = false;
			bool ti1z = false;
			bool ti2z = false;
			bool ti3z = false;
			bool ti4z = false;

			if (rz[0] && rz[2]) {
				tr1z = true;
				tr3z = true;
			} else if (!rz[0] && rz[2]) {
			} else if (rz[0] && !rz[2]) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (rz[1] && rz[3]) {
				tr2z = true;
				tr4z = true;
			} else if (!rz[1] && rz[3]) {
			} else if (rz[1] && !rz[3]) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (iz[0] && iz[2]) {
				ti1z = true;
				ti3z = true;
			} else if (!iz[0] && iz[2]) {
			} else if (iz[0] && !iz[2]) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (iz[1] && iz[3]) {
				ti2z = true;
				ti4z = true;
			} else if (!iz[1] && iz[3]) {
			} else if (iz[1] && !iz[3]) {
				cnt++;
			} else {
				cnt += 2;
			}
			std::unordered_set<int> used;
			if (tr1z && tr2z) {
			} else if (!tr1z && tr2z) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
				}
				if (used.find(k2) == used.end()) {
					used.insert(k2);
				}
			} else if (tr1z && !tr2z) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
				}
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					cnt++;
				}
			} else {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					cnt++;
				}
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					cnt++;
				}
			}

			if (tr3z && ti4z) {
			} else if (!tr3z && ti4z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
				}
			} else if (tr3z && !ti4z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					cnt++;
				}
			} else {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					cnt++;
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					cnt++;
				}
			}
			used.clear();
			if (ti1z && ti2z) {
			} else if (!ti1z && ti2z) {
				if (used.find(k0) == used.end() && k0 != N - k0 && k0 != 0) {
					used.insert(k0);
				}
				if (used.find(k2) == used.end() && k2 != N - k2 && k2 != 0) {
					cnt++;
					used.insert(k2);
				}
			} else if (ti1z && !ti2z) {
				if (used.find(k0) == used.end() && k0 != N - k0 && k0 != 0) {
					used.insert(k0);
				}
				if (used.find(k2) == used.end() && k2 != N - k2 && k2 != 0) {
					used.insert(k2);
				}
			} else {
				if (used.find(k0) == used.end() && k0 != N - k0 && k0 != 0) {
					cnt++;
					used.insert(k0);
				}
				if (used.find(k2) == used.end() && k2 != N - k2 && k2 != 0) {
					cnt += 2;
					used.insert(k2);
				}
			}
			if (ti3z && tr4z) {
			} else if (!ti3z && tr4z) {
				if (used.find(k1) == used.end() && k1 != N - k1 && k1 != 0) {
					used.insert(k1);
				}
				if (used.find(k3) == used.end() && k3 != N - k3 && k3 != 0) {
					cnt++;
					used.insert(k3);
				}
			} else if (ti3z && !tr4z) {
				if (used.find(k1) == used.end() && k1 != N - k1 && k1 != 0) {
					cnt++;
					used.insert(k1);
				}
				if (used.find(k3) == used.end() && k3 != N - k3 && k3 != 0) {
					cnt++;
					used.insert(k3);
				}
			} else {
				if (used.find(k1) == used.end() && k1 != N - k1 && k1 != 0) {
					cnt++;
					used.insert(k1);
				}
				if (used.find(k3) == used.end() && k3 != N - k3 && k3 != 0) {
					cnt += 2;
					used.insert(k3);
				}
			}
		}
			break;
		case 5: {
			// Radix -- 5
			bool tr1z = false;
			bool tr2z = false;
			bool tr3z = false;
			bool tr4z = false;
			bool tr5z = false;
			bool tr6z = false;
			bool tr7z = false;
			bool tr8z = false;
			bool tr9z = false;
			bool tr10z = false;
			bool tr11z = false;
			bool ti1z = false;
			bool ti2z = false;
			bool ti3z = false;
			bool ti4z = false;
			bool ti5z = false;
			bool ti6z = false;
			bool ti7z = false;
			bool ti8z = false;
			bool ti9z = false;
			bool ti10z = false;
			bool ti11z = false;
			if (rz[1] && rz[4]) {
				tr1z = tr3z = true;
			} else if (!rz[1] && rz[4]) {
			} else if (rz[1] && !rz[4]) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (rz[2] && rz[3]) {
				tr2z = tr4z = true;
			} else if (!rz[2] && rz[3]) {
			} else if (rz[2] && !rz[3]) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (tr1z && tr2z) {
				tr5z = true;
				tr6z = true;
			} else if (!tr1z && tr2z) {
				cnt++;
			} else if (tr2z && !tr1z) {
				cnt++;
			} else {
				cnt += 3;
			}
			if (rz[0] && tr5z) {
				tr7z = true;
			} else if (!rz[0] && tr5z) {
			} else if (rz[0] && !tr5z) {
				cnt++;
			} else {
				cnt++;
			}
			if (tr6z && tr7z) {
				tr8z = true;
				tr9z = true;
			} else if (tr6z && !tr7z) {
			} else if (!tr6z && tr7z) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (tr3z && tr4z) {
				tr10z = tr11z = true;
			} else if (!tr3z && tr4z) {
				cnt += 2;
			} else if (tr3z && !tr4z) {
				cnt += 2;
			} else {
				cnt += 4;
			}

			if (iz[1] && iz[4]) {
				ti1z = ti3z = true;
			} else if (!iz[1] && iz[4]) {
			} else if (iz[1] && !iz[4]) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (iz[2] && iz[3]) {
				ti2z = ti4z = true;
			} else if (!iz[2] && iz[3]) {
			} else if (iz[2] && !iz[3]) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (ti1z && ti2z) {
				ti5z = true;
				ti6z = true;
			} else if (!ti1z && ti2z) {
				cnt++;
			} else if (ti2z && !ti1z) {
				cnt++;
			} else {
				cnt += 3;
			}
			if (iz[0] && ti5z) {
				ti7z = true;
			} else if (!iz[0] && ti5z) {
			} else if (iz[0] && !ti5z) {
				cnt++;
			} else {
				cnt++;
			}
			if (ti6z && ti7z) {
				ti8z = true;
				ti9z = true;
			} else if (ti6z && !ti7z) {
			} else if (!ti6z && ti7z) {
				cnt++;
			} else {
				cnt += 2;
			}
			if (ti3z && ti4z) {
				ti10z = ti11z = true;
			} else if (!ti3z && ti4z) {
				cnt += 2;
			} else if (ti3z && !ti4z) {
				cnt += 2;
			} else {
				cnt += 4;
			}
			int k0 = k;
			int k1 = k + N / 5;
			int k2 = k + 2 * N / 5;
			int k3 = N - (k + 3 * N / 5);
			int k4 = N - (k + 4 * N / 5);
			if (k > N / 2) {
				k2 = N - k2;
			}
			std::unordered_set<int> used;
			if (used.find(k0) == used.end()) {
				used.insert(k0);
				if (rz[0] && tr5z) {
				} else if (!rz[0] && tr5z) {
				} else if (rz[0] && !tr5z) {
				} else {
					cnt++;
				}
			}
			if (tr8z && ti10z) {
			} else if (!tr8z && ti10z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
				}
			} else if (tr8z && !ti10z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
					cnt++;
				}
			} else {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					cnt++;
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
					cnt++;
				}
			}
			if (tr9z && ti11z) {
			} else if (!tr9z && ti11z) {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
				}
			} else if (tr9z && !ti11z) {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					cnt++;
				}
			} else {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					cnt++;
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					cnt++;
				}
			}

			used.clear();
			if (used.find(k0) == used.end() && k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
				used.insert(k0);
				if (iz[0] && ti5z) {
				} else if (!iz[0] && ti5z) {
				} else if (iz[0] && !ti5z) {
				} else {
					cnt++;
				}
			}
			if (ti8z && tr10z) {
			} else if (!ti8z && tr10z) {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
					cnt++;
				}
			} else if (ti8z && !tr10z) {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
					cnt += 2;
				}
			} else {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
					cnt += 3;
				}
			}
			if (ti9z && tr11z) {
			} else if (!ti9z && tr11z) {
				if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					used.insert(k2);
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					cnt++;
				}
			} else if (ti9z && !tr11z) {
				if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					used.insert(k2);
					cnt++;
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					cnt++;
				}
			} else {
				if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
					used.insert(k2);
					cnt++;
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					cnt += 2;
				}
			}
		}
			break;
		case 6: {
			bool tr1z = false;
			bool tr2z = false;
			bool tr3z = false;
			bool tr4z = false;
			bool tr5z = false;
			bool tr6z = false;
			bool tr7z = false;
			bool tr8z = false;
			bool tr9z = false;
			bool tr10z = false;
			bool tr11z = false;
			bool tr12z = false;
			bool ti1z = false;
			bool ti2z = false;
			bool ti3z = false;
			bool ti4z = false;
			bool ti5z = false;
			bool ti6z = false;
			bool ti7z = false;
			bool ti8z = false;
			bool ti9z = false;
			bool ti10z = false;
			bool ti11z = false;
			bool ti12z = false;

			if (rz[2] && rz[4]) {
				tr1z = true;
				tr3z = true;
			} else if (!rz[2] && rz[4]) {
				cnt++;
			} else if (rz[2] && !rz[4]) {
				cnt++;
			} else {
				cnt += 3;
			}
			if (tr1z && rz[0]) {
				tr2z = true;
			} else if (!tr1z && rz[0]) {
				cnt++;
			} else if (tr1z && !rz[0]) {
			} else {
				cnt++;
			}
			if (rz[1] && rz[5]) {
				tr4z = true;
				tr6z = true;
			} else if (!rz[1] && rz[5]) {
				cnt++;
			} else if (rz[1] && !rz[5]) {
				cnt++;
			} else {
				cnt += 3;
			}
			if (tr4z && rz[3]) {
				tr5z = true;
			} else if (!tr4z && rz[3]) {
				cnt++;
			} else if (tr4z && !rz[3]) {
			} else {
				cnt++;
			}

			if (iz[2] && iz[4]) {
				ti1z = true;
				ti3z = true;
			} else if (!iz[2] && iz[4]) {
				cnt++;
			} else if (iz[2] && !iz[4]) {
				cnt++;
			} else {
				cnt += 3;
			}
			if (ti1z && iz[0]) {
				ti2z = true;
			} else if (!ti1z && iz[0]) {
				cnt++;
			} else if (ti1z && !iz[0]) {
			} else {
				cnt++;
			}
			if (iz[1] && iz[5]) {
				ti4z = true;
				ti6z = true;
			} else if (!iz[1] && iz[5]) {
				cnt++;
			} else if (iz[1] && !iz[5]) {
				cnt++;
			} else {
				cnt += 3;
			}
			if (ti4z && iz[3]) {
				ti5z = true;
			} else if (!ti4z && iz[3]) {
				cnt++;
			} else if (tr4z && !rz[3]) {
			} else {
				cnt++;
			}

			if (rz[0] && tr1z) {
				tr7z = true;
			} else if (!rz[0] && tr1z) {
			} else if (rz[0] && !tr1z) {
			} else {
				cnt++;
			}

			if (iz[0] && ti1z) {
				ti7z = true;
			} else if (!iz[0] && ti1z) {
			} else if (iz[0] && !ti1z) {
			} else {
				cnt++;
			}

			if (tr2z && ti3z) {
				tr8z = true;
				tr9z = true;
			} else if (!tr2z && ti3z) {
			} else if (tr2z && !ti3z) {
				cnt++;
			} else {
				cnt += 2;
			}

			if (ti2z && tr3z) {
				ti2z = true;
				tr3z = true;
			} else if (!ti2z && tr3z) {
			} else if (ti2z && !tr3z) {
				cnt++;
			} else {
				cnt += 2;
			}

			if (rz[3] && tr4z) {
				tr10z = true;
			} else if (rz[3] && !tr4z) {
			} else if (!rz[3] && tr4z) {
			} else {
				cnt++;
			}

			if (iz[3] && ti4z) {
				ti4z = true;
			} else if (iz[3] && !ti4z) {
			} else if (!iz[3] && ti4z) {
			} else {
				cnt++;
			}

			if (tr5z && ti6z) {
				tr5z = true;
				ti6z = true;
			} else if (!tr5z && ti6z) {
			} else if (tr5z && !ti6z) {
				cnt++;
			} else {
				cnt += 2;
			}

			if (ti5z && tr6z) {
				ti5z = true;
				tr6z = true;
			} else if (!ti5z && tr6z) {
			} else if (ti5z && !tr6z) {
				cnt++;
			} else {
				cnt += 2;
			}

			int k0 = k + 0 * N / 6;
			int k1 = k + 1 * N / 6;
			int k2 = k + 2 * N / 6;
			int k3 = N - (k + 3 * N / 6);
			int k4 = N - (k + 4 * N / 6);
			int k5 = N - (k + 5 * N / 6);
			std::unordered_set<int> used;
			if (tr7z && tr10z) {
			} else if (!tr7z && tr10z) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
				}
			} else if (tr7z && !tr10z) {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					cnt++;
				}
			} else {
				if (used.find(k0) == used.end()) {
					used.insert(k0);
					cnt++;
				}
				if (used.find(k3) == used.end()) {
					used.insert(k3);
					cnt++;
				}
			}
			if (tr8z && tr11z) {

			} else if (!tr8z && tr11z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
				}
			} else if (tr8z && !tr11z) {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					cnt++;
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
				}
			} else {
				if (used.find(k1) == used.end()) {
					used.insert(k1);
					cnt++;
				}
				if (used.find(k4) == used.end()) {
					used.insert(k4);
					cnt++;
				}
			}
			if (tr9z && tr12z) {
			} else if (tr9z && !tr12z) {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
				}
				if (used.find(k5) == used.end()) {
					used.insert(k5);
					cnt++;
				}
			} else if (!tr9z && tr12z) {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
				}
				if (used.find(k5) == used.end()) {
					used.insert(k5);
				}
			} else {
				if (used.find(k2) == used.end()) {
					used.insert(k2);
					cnt++;
				}
				if (used.find(k5) == used.end()) {
					used.insert(k5);
					cnt++;
				}
			}
			used.clear();

			if (ti7z && ti10z) {
			} else if (!ti7z && ti10z) {
				if (used.find(k0) == used.end() && k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
					used.insert(k0);
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					cnt++;
				}
			} else if (ti7z && !ti10z) {
				if (used.find(k0) == used.end() && k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
					used.insert(k0);
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
				}
			} else {
				if (used.find(k0) == used.end() && k0 != 0 && !(k0 == N / 2 && N % 2 == 0)) {
					used.insert(k0);
					cnt++;
				}
				if (used.find(k3) == used.end() && k3 != 0 && !(k3 == N / 2 && N % 2 == 0)) {
					used.insert(k3);
					cnt++;
				}
			}

			if (ti8z && ti11z) {

			} else if (!ti8z && ti11z) {
				if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
					used.insert(k1);
				}
				if (used.find(k4) == used.end() && k4 != 0 && !(k4 == N / 2 && N % 2 == 0)) {
					used.insert(k4);
					cnt++;

				} else if (ti8z && !ti11z) {
					if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
						used.insert(k1);
						cnt++;
					}
					if (used.find(k4) == used.end() && k4 != 0 && !(k4 == N / 2 && N % 2 == 0)) {
						used.insert(k4);
						cnt++;
					}
				} else {
					if (used.find(k1) == used.end() && k1 != 0 && !(k1 == N / 2 && N % 2 == 0)) {
						used.insert(k1);
						cnt++;
					}
					if (used.find(k4) == used.end() && k4 != 0 && !(k4 == N / 2 && N % 2 == 0)) {
						used.insert(k4);
						cnt += 2;
					}
				}
				if (ti9z && ti12z) {
				} else if (!ti9z && ti12z) {
					if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
						used.insert(k2);
					}
					if (used.find(k5) == used.end() && k5 != 0 && !(k5 == N / 2 && N % 2 == 0)) {
						used.insert(k5);
						cnt++;
					}
				} else if (ti9z && !ti12z) {
					if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
						used.insert(k2);
					}
					if (used.find(k5) == used.end() && k5 != 0 && !(k5 == N / 2 && N % 2 == 0)) {
						used.insert(k5);
					}
				} else {
					if (used.find(k2) == used.end() && k2 != 0 && !(k2 == N / 2 && N % 2 == 0)) {
						used.insert(k2);
						cnt++;
					}
					if (used.find(k5) == used.end() && k5 != 0 && !(k5 == N / 2 && N % 2 == 0)) {
						used.insert(k5);
						cnt++;
					}
				}
			}
		}
			break;
		default: {
			std::vector<bool> txpz((r - 1) / 2 + 1, false);
			std::vector<bool> txmz((r - 1) / 2 + 1, false);
			std::vector<bool> typz((r - 1) / 2 + 1, false);
			std::vector<bool> tymz((r - 1) / 2 + 1, false);
			for (int j = 1; j <= (r - 1) / 2; j++) {
				if (rz[j] && rz[r - j]) {
					txpz[j] = true;
				} else if (!rz[j] && rz[r - j]) {
				} else if (rz[j] && !rz[r - j]) {
				} else {
					cnt++;
				}
			}
			for (int j = 1; j <= (r - 1) / 2; j++) {
				if (rz[j] && rz[r - j]) {
					txmz[j] = true;
				} else if (!rz[j] && rz[r - j]) {
				} else if (rz[j] && !rz[r - j]) {
					cnt++;
				} else {
					cnt++;
				}
			}
			for (int j = 1; j <= (r - 1) / 2; j++) {
				if (iz[j] && iz[r - j]) {
					typz[j] = true;
				} else if (!iz[j] && iz[r - j]) {
				} else if (iz[j] && !iz[r - j]) {
				} else {
					cnt++;
				}
			}
			for (int j = 1; j <= (r - 1) / 2; j++) {
				if (iz[j] && iz[r - j]) {
					tymz[j] = true;
				} else if (!iz[j] && iz[r - j]) {
				} else if (iz[j] && !iz[r - j]) {
					cnt++;
				} else {
					cnt++;
				}
			}
			std::vector<bool> apz((r - 1) / 2 + 1, true);
			std::vector<bool> amz((r - 1) / 2 + 1, true);
			std::vector<bool> bpz((r - 1) / 2 + 1, true);
			std::vector<bool> bmz((r - 1) / 2 + 1, true);
			for (int i = 1; i <= (r - 1) / 2; i++) {
				bool apfirst = true;
				bool amfirst = true;
				bool bpfirst = true;
				bool bmfirst = true;

				for (int j = 1; j <= (r - 1) / 2; j++) {
					if (!txpz[j]) {
						apz[i] = false;
						if (apfirst) {
							apfirst = false;
							cnt++;
						} else {
							cnt++;
						}
					}
					if (!typz[j]) {
						bpz[i] = false;
						if (bpfirst) {
							bpfirst = false;
							cnt++;
						} else {
							cnt++;
						}
					}
					if (!tymz[j]) {
						amz[i] = false;
						if (amfirst) {
							cnt++;
							amfirst = false;
						} else {
							cnt++;
						}
					}
					if (!txmz[j]) {
						bmz[i] = false;
						if (bmfirst) {
							cnt++;
							bmfirst = false;
						} else {
							cnt++;
						}
					}
				}
				if (!rz[0]) {
					apz[i] = false;
					if (apfirst) {
						apfirst = false;
					} else {
						cnt++;
					}
				}
				if (!iz[0]) {
					bpz[i] = false;
					if (bpfirst) {
						bpfirst = false;
					} else {
						cnt++;
					}
				}
			}
			std::unordered_set<int> used;
			if (used.find(k) == used.end()) {
				used.insert(k);

				for (int i = 0; i < r; i++) {
					if (!rz[i]) {
						cnt++;
					}
				}

			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				const int ka = k + i * N / r;
				const int rk0 = index_real(o, ka, 0, N);
				if (used.find(ka) == used.end()) {
					used.insert(ka);
					if (apz[i] && amz[i]) {
					} else if (apz[i] && !amz[i]) {
					} else if (!apz[i] && amz[i]) {
					} else {
						cnt++;
					}
				}
			}
			if (!(k == N / r / 2 && N % 2 == 0)) {
				for (int i = 1; i <= (r - 1) / 2; i++) {
					const int kb = N - (k + (r - i) * N / r);
					const int rk1 = index_real(o, kb, 0, N);
					if (used.find(kb) == used.end()) {
						used.insert(kb);
						if (apz[i] && amz[i]) {
						} else if (apz[i] && !amz[i]) {
							cnt++;
						} else if (!apz[i] && amz[i]) {
						} else {
							cnt++;
						}
					}
				}
			}
			used.clear();
			if (used.find(k) == used.end()) {
				if (k != 0 && !(k == N / 2 && N % 2 == 0)) {
					used.insert(k);
					for (int i = 0; i < r; i++) {
						if (!iz[i]) {
							cnt++;
						}
					}
				}
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				const int ka = k + i * N / r;
				const int ik0 = index_real(o, ka, 1, N);
				if (used.find(ka) == used.end() && ka != 0 && !(ka == N / 2 && N % 2 == 0)) {
					used.insert(ka);
					if (bpz[i] && bmz[i]) {
					} else if (bpz[i] && !bmz[i]) {
						cnt++;
					} else if (!bpz[i] && bmz[i]) {
					} else {
						cnt++;
					}
				}
			}
			for (int i = 1; i <= (r - 1) / 2; i++) {
				const int kb = N - (k + (r - i) * N / r);
				const int ik1 = index_real(o, kb, 1, N);
				if (used.find(kb) == used.end() && kb != 0 && !(kb == N / 2 && N % 2 == 0)) {
					used.insert(kb);
					if (bpz[i] && bmz[i]) {
					} else if (bpz[i] && !bmz[i]) {
						cnt++;
					} else if (!bpz[i] && bmz[i]) {
						cnt++;
					} else {
						cnt += 2;
					}
				}
			}
		}
		}
	}

	return cnt;
}
