#include "fftgen.hpp"

void print_short_fft(int N, std::vector<std::string> in_, std::vector<std::string> out_, bool real = false) {
	print("/* begin short transform - %s - length %i */\n", real ? "real" : "complex", N);
	std::vector<const char*> in(in_.size());
	std::vector<const char*> out(out_.size());
	for (int i = 0; i < in.size(); i++) {
		in[i] = in_[i].c_str();
		out[i] = out_[i].c_str();
	}

	auto pfac = prime_factorization(N);
	if (pfac.size() == 1) {
		std::vector<int> ri(N);
		std::vector<int> ii(N);
		if (real) {
			for (int n = 0; n < N; n++) {
				ri[n] = std::min(n, N - n);
				ii[n] = std::max(n, N - n);
			}
		} else {
			for (int n = 0; n < N; n++) {
				ri[n] = 2 * n;
				ii[n] = 2 * n + 1;
			}
		}
		std::vector<std::string> vars;
		if (N > 2) {
			for (int n = 1; n <= (N - 1) / 2; n++) {
				vars.push_back(std::string("ap") + std::to_string(n));
				vars.push_back(std::string("bm") + std::to_string(n));
				if (!real) {
					vars.push_back(std::string("bp") + std::to_string(n));
					vars.push_back(std::string("am") + std::to_string(n));
				}
			}
		}
		vars.push_back("sr0");
		vars.push_back("si0");
		if (N % 2 == 0) {
			vars.push_back("pr0");
			vars.push_back("pi0");
			if (N > 4) {
				vars.push_back("er0");
				vars.push_back("ei0");
				vars.push_back("or0");
				vars.push_back("oi0");
			}
		}
		if (vars.size()) {
			int cnt = 1;
			int comma = true;
			print("double %s", vars.front().c_str());
			for (int i = 1; i < vars.size(); i++) {
				if (comma) {
					print_notab(", ");
				}
				print_notab("%s", vars[i].c_str());
				cnt++;
				if (cnt % 8 == 0) {
					if (i == vars.size() - 1) {
					} else {
						print_notab(",\n");
						print("       ");
					}
					comma = false;
				} else {
					comma = true;
				}
			}
			print_notab(";\n");
		}
		if (!real) {
			for (int j = 1; j <= (N - 1) / 2; j++) {
				print("const auto tp%i = %s + %s;\n", j, in[ri[j]], in[ri[N - j]]);
			}
			for (int j = 1; j <= (N - 1) / 2; j++) {
				print("const auto tm%i = %s - %s;\n", j, in[ri[j]], in[ri[N - j]]);
			}
			for (int j = 1; j <= (N - 1) / 2; j++) {
				print("const auto up%i = %s + %s;\n", j, in[ii[j]], in[ii[N - j]]);
			}
			for (int j = 1; j <= (N - 1) / 2; j++) {
				print("const auto um%i = %s - %s;\n", j, in[ii[j]], in[ii[N - j]]);
			}
		} else {
			for (int j = 1; j <= (N - 1) / 2; j++) {
				print("const auto tp%i = %s + %s;\n", j, in[ri[j]], in[ii[N - j]]);
			}
			for (int j = 1; j <= (N - 1) / 2; j++) {
				print("const auto tm%i = %s - %s;\n", j, in[ri[j]], in[ii[N - j]]);
			}
		}
		for (int i = 1; i <= (N - 1) / 2; i++) {
			print("ap%i = %s;\n", i, in[ri[0]]);
			if (!real) {
				print("bp%i = %s;\n", i, in[ii[0]]);
			}
			for (int j = 1; j <= (N - 1) / 2; j++) {
				int k = (j * i) % N;
				if (k == 0) {
					print("ap%i += tp%i;\n", i, j);
					real ? void() : print("bp%i += up%i;\n", i, j);
				} else if (k == N / 2 && N % 2 == 0) {
					print("ap%i -= tp%i;\n", i, j);
					real ? void() : print("bp%i -= up%i;\n", i, j);
				} else if (k == N / 4 && N % 4 == 0) {
				} else if (k == 3 * N / 4 && N % 4 == 0) {
				} else {
					print("ap%i = std::fma(tp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / N), i);
					real ? void() : print("bp%i = std::fma(up%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / N), i);
				}
			}
			if (N % 2 == 0) {
				print("ap%i %c= %s;\n", i, i % 2 == 0 ? '+' : '-', in[ri[N / 2]]);
				real ? void() : print("bp%i %c= %s;\n", i, i % 2 == 0 ? '+' : '-', in[ii[N / 2]]);
			}
		}
		for (int i = 1; i <= (N - 1) / 2; i++) {
			bool init = false;
			std::vector<bool> init_ele((N - 1) / 2 + 1, false);
			for (int j = 1; j <= (N - 1) / 2; j++) {
				int k = (j * i) % N;
				if (k == N / 4 && N % 4 == 0) {
					real ? void() : print("am%i = um%i;\n", i, j);
					print("bm%i = tm%i;\n", i, j);
					init = true;
					init_ele[j] = true;
					break;
				}
			}
			if (!init) {
				for (int j = 1; j <= (N - 1) / 2; j++) {
					int k = (j * i) % N;
					if (k == 0) {
					} else if (k == N / 2 && N % 2 == 0) {
					} else if (k == N / 4 && N % 4 == 0) {
					} else if (k == 3 * N / 4 && N % 4 == 0) {
						real ? void() : print("am%i = -um%i;\n", i, j);
						print("bm%i = -tm%i;\n", i, j);
						init = true;
						init_ele[j] = true;
						break;
					} else {
						real ? void() : print("am%i = um%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N));
						print("bm%i = tm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N));
						init = true;
						init_ele[j] = true;
						break;
					}
				}
			}
			for (int j = 1; j <= (N - 1) / 2; j++) {
				if (!init_ele[j]) {
					int k = (j * i) % N;
					if (k == 0) {
					} else if (k == N / 2 && N % 2 == 0) {
					} else if (k == N / 4 && N % 4 == 0) {
						real ? void() : print("am%i += um%i;\n", i, j);
						print("bm%i += tm%i;\n", i, j);
					} else if (k == 3 * N / 4 && N % 4 == 0) {
						real ? void() : print("am%i -= um%i;\n", i, j);
						print("bm%i -= tm%i;\n", i, j);
					} else {
						real ? void() : print("am%i = std::fma(um%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / N), i);
						print("bm%i = std::fma(tm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / N), i);
					}
				}
			}
		}
		if (N % 2 == 0) {
			if (N % 4 == 0) {
				print("pr0 = sr0 = %s + %s;\n", in[ri[0]], in[ri[N / 2]]);
				if (!real) {
					print("pi0 = si0 = %s + %s;\n", in[ii[0]], in[ii[N / 2]]);
				}
			} else {
				print("sr0 = %s + %s;\n", in[ri[0]], in[ri[N / 2]]);
				print("pr0 = %s - %s;\n", in[ri[0]], in[ri[N / 2]]);
				if (!real) {
					print("si0 = %s + %s;\n", in[ii[0]], in[ii[N / 2]]);
					print("pi0 = %s - %s;\n", in[ii[0]], in[ii[N / 2]]);
				}
			}
		} else {
			print("sr0 = %s;\n", in[ri[0]]);
			if (!real) {
				print("si0 = %s;\n", in[ii[0]]);
			}
		}
		if (N % 2 == 0 && N > 4) {
			print("er0 = tp2;\n");
			print("or0 = tp1;\n");
			if (!real) {
				print("ei0 = up2;\n");
				print("oi0 = up1;\n");
			}
			for (int i = 4; i <= (N - 1) / 2; i += 2) {
				print("er0 += tp%i;\n", i);
			}
			for (int i = 3; i <= (N - 1) / 2; i += 2) {
				print("or0 += tp%i;\n", i);
			}
			print("sr0 += er0 + or0;\n");
			print("pr0 += er0 - or0;\n");
			if (!real) {
				for (int i = 4; i <= (N - 1) / 2; i += 2) {
					print("ei0 += up%i;\n", i);
				}
				for (int i = 3; i <= (N - 1) / 2; i += 2) {
					print("oi0 += up%i;\n", i);
				}
				print("si0 += ei0 + oi0;\n");
				print("pi0 += ei0 - oi0;\n");
			}
		} else {
			for (int i = 1; i <= (N - 1) / 2; i++) {
				print("sr0 += tp%i;\n", i);
			}
			if (!real) {
				for (int i = 1; i <= (N - 1) / 2; i++) {
					print("si0 += up%i;\n", i);
				}
			}
			if (N % 2 == 0) {
				for (int i = 1; i <= (N - 1) / 2; i++) {
					print("pr0 %c= tp%i;\n", i % 2 == 0 ? '+' : '-', i);
				}
				if (!real) {
					for (int i = 1; i <= (N - 1) / 2; i++) {
						print("pi0 %c= up%i;\n", i % 2 == 0 ? '+' : '-', i);
					}
				}
			}
		}
		if (!real) {
			print("%s = sr0;\n", out[ri[0]]);
			print("%s = si0;\n", out[ii[0]]);
			if (N % 2 == 0) {
				print("%s = pr0;\n", out[ri[N / 2]]);
				print("%s = pi0;\n", out[ii[N / 2]]);
			}
			for (int i = 1; i <= (N - 1) / 2; i++) {
				print("%s = ap%i + am%i;\n", out[ri[i]], i, i);
				print("%s = bp%i - bm%i;\n", out[ii[i]], i, i);
				print("%s = ap%i - am%i;\n", out[ri[N - i]], i, i);
				print("%s = bp%i + bm%i;\n", out[ii[N - i]], i, i);
			}
		} else {
			print("%s = sr0;\n", out[ri[0]]);
			if (N % 2 == 0) {
				print("%s = pr0;\n", out[ri[N / 2]]);
			}
			for (int i = 1; i <= (N - 1) / 2; i++) {
				print("%s = ap%i;\n", out[ri[i]], i, i);
				print("%s = -bm%i;\n", out[ii[i]], i, i);
			}
		}
	} else if (real) {
		int N2 = pow(pfac[0].first, pfac[0].second);
		int N1 = N / N2;
		std::vector<std::vector<std::string>> in1(N2, std::vector<std::string>(N1));
		std::vector<std::vector<std::string>> out1(N2, std::vector<std::string>(N1));
		std::vector<std::vector<std::string>> out2(N1 / 2 + 1, std::vector<std::string>(2 * N2));
		std::vector<std::vector<std::string>> in2(N1 / 2 + 1, std::vector<std::string>(2 * N2));
		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				const int nn = ((n1 * N2 + n2 * N1) % N);
				in1[n2][n1] = in[nn];
			}
		}
		print("double zr0, ");
		for (int n = 1; n < N - 1; n++) {
			print_notab("zr%i, ", n);
		}
		print_notab("zr%i;\n", N - 1);
		print("double ");
		for (int n = 1; n < N - 2; n++) {
			print_notab("zi%i, ", n);
		}
		print_notab("zi%i;\n", N - 2);
		for (int n2 = 0; n2 < N2; n2++) {
			out1[n2][0] = std::string("zr") + std::to_string(n2 * N1);
			for (int n1 = 1; n1 < N1 / 2 + 1; n1++) {
				const int nn = ((n1 * N2 + n2 * N1) % N);
				out1[n2][n1] = std::string("zr") + std::to_string(n2 * N1 + n1);
				if (n1 != N1 - n1) {
					out1[n2][N1 - n1] = std::string("zi") + std::to_string(n2 * N1 + n1);
				}
			}
		}
		for (int n2 = 0; n2 < N2; n2++) {
			in2[0][n2] = out1[n2][0];
			for (int n1 = 1; n1 < N1 / 2 + 1; n1++) {
				in2[n1][2 * n2] = out1[n2][n1];
				in2[n1][2 * n2 + 1] = out1[n2][N1 - n1];
			}
		}
		std::vector<std::string> negs;
		for (int n = 0; n < N; n++) {
			int n1 = n % N1;
			int n2 = n % N2;
			if (n1 == 0) {
				if (n2 > N2 / 2) {
					out2[n1][n2] = out[std::max(n, N - n)];
					negs.push_back(out[std::max(n, N - n)]);
				} else {
					out2[n1][n2] = out[std::min(n, N - n)];
				}
			} else if (n1 < N1 / 2 + 1) {
				out2[n1][2 * n2] = out[std::min(n, N - n)];
				if (n != 0) {
					out2[n1][2 * n2 + 1] = out[std::max(n, N - n)];
					if (n > N - n) {
						negs.push_back(out[n]);
					}
				}
			}

		}
		for (int n2 = 0; n2 < N2; n2++) {
			print("{\n");
			indent();
			print_real_short_fft(N1, in1[n2], out1[n2]);
			deindent();
			print("}\n");
		}
		{
			int n1 = 0;
			print("{\n");
			indent();
			print_real_short_fft(N2, in2[n1], out2[n1]);
			deindent();
			print("}\n");
		}
		for (int n1 = 1; n1 < N1 / 2 + 1; n1++) {
			print("{\n");
			indent();
			print_complex_short_fft(N2, in2[n1], out2[n1]);
			deindent();
			print("}\n");
		}
		for (int l = 0; l < negs.size(); l++) {
			print("%s = -%s;\n", negs[l].c_str(), negs[l].c_str());
		}
	} else {
		int N1 = pow(pfac[0].first, pfac[0].second);
		int N2 = N / N1;
		if (N == 12) {
			printf("12 =  %i %i\n", N1, N2);
		}

		print("double qr0, ");
		for (int n = 1; n < N - 1; n++) {
			print_notab("qr%i, ", n);
		}
		print_notab("qr%i;\n", N - 1);
		print("double qi0, ");
		for (int n = 1; n < N - 1; n++) {
			print_notab("qi%i, ", n);
		}
		print_notab("qi%i;\n", N - 1);
		std::vector<std::vector<std::string>> out1(N2, std::vector<std::string>(2 * N1));
		std::vector<std::vector<std::string>> in1(N2, std::vector<std::string>(2 * N1));
		std::vector<std::vector<std::string>> out2(N1, std::vector<std::string>(2 * N2));
		std::vector<std::vector<std::string>> in2(N1, std::vector<std::string>(2 * N2));
		for (int n1 = 0; n1 < N1; n1++) {
			for (int n2 = 0; n2 < N2; n2++) {
				const int nn = ((n1 * N2 + n2 * N1) % N);
				in1[n2][2 * n1] = in[2 * nn];
				in1[n2][2 * n1 + 1] = in[2 * nn + 1];
				out1[n2][2 * n1] = std::string("qr") + std::to_string(n1 + n2 * N1);
				out1[n2][2 * n1 + 1] = std::string("qi") + std::to_string(n1 + n2 * N1);
				in2[n1][2 * n2] = std::string("qr") + std::to_string(n1 + n2 * N1);
				in2[n1][2 * n2 + 1] = std::string("qi") + std::to_string(n1 + n2 * N1);
			}
		}
		for (int n = 0; n < N; n++) {
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
	print("/* end short transform - %s - length %i */\n", real ? "real" : "complex", N);

}

void print_complex_short_fft(int N, std::vector<std::string> in_, std::vector<std::string> out_) {
	print_short_fft(N, in_, out_);
}

void print_real_short_fft(int N, std::vector<std::string> in_, std::vector<std::string> out_) {
	print_short_fft(N, in_, out_, true);
}

/*
 void print_real_short_fft(int r, std::vector<std::string> in_, std::vector<std::string> out_) {
 std::vector<const char*> in(in_.size());
 std::vector<const char*> out(out_.size());
 for (int i = 0; i < in.size(); i++) {
 in[i] = in_[i].c_str();
 out[i] = out_[i].c_str();
 }
 switch (r) {
 case 2:
 print("const auto tr = %s;\n", in[0]);
 print("%s = tr + %s;\n", out[0], in[1]);
 print("%s = tr - %s;\n", out[1], in[1]);
 break;
 case 3:
 print("const auto tr1 = %s + %s;\n", in[1], in[2]);
 print("const auto tr2 = std::fma(tr1, -0.5, %s);\n", in[0]);
 print("const auto tr3 = (%24.17e) * (%s - %s);\n", sqrt(3) * 0.5, in[1], in[2]);
 print("%s = %s + tr1;\n", out[0], in[0]);
 print("%s = tr2;\n", out[1]);
 print("%s = -tr3;\n", out[2]);
 break;
 case 4:
 print("const auto tr1 = %s + %s;\n", in[0], in[2]);
 print("const auto tr3 = %s - %s;\n", in[0], in[2]);
 print("const auto tr2 = %s + %s;\n", in[1], in[3]);
 print("const auto tr4 = %s - %s;\n", in[1], in[3]);
 print("%s = tr1 + tr2;\n", out[0]);
 print("%s = tr3;\n", out[1]);
 print("%s = tr1 - tr2;\n", out[2]);
 print("%s = -tr4;\n", out[3]);
 break;
 case 5:
 print("const auto tr1 = %s + %s;\n", in[1], in[4]);
 print("const auto tr2 = %s + %s;\n", in[2], in[3]);
 print("const auto tr3 = %s - %s;\n", in[1], in[4]);
 print("const auto tr4 = %s - %s;\n", in[2], in[3]);
 print("const auto tr5 = tr1 + tr2;\n");
 print("const auto tr6 = (%24.17e) * (tr1 - tr2);\n", sqrt(5) * 0.25);
 print("const auto tr7 = std::fma(tr5, -0.25, %s);\n", in[0]);
 print("const auto tr8 = tr7 + tr6;\n");
 print("const auto tr9 = tr7 - tr6;\n");
 print("const auto tr10 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
 print("const auto tr11 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
 print("%s = %s + tr5;\n", out[0], in[0]);
 print("%s = tr8;\n", out[1]);
 print("%s = -tr10;\n", out[4]);
 print("%s = tr9;\n", out[2]);
 print("%s = -tr11;\n", out[3]);
 break;
 case 6:
 print("const auto tr1 = %s + %s;\n", in[2], in[4]);
 print("const auto tr2 = std::fma(tr1, -0.5, %s);\n", in[0]);
 print("const auto tr3 = (%24.17e) * (%s - %s);\n", sin(M_PI / 3.0), in[2], in[4]);
 print("const auto tr4 = %s + %s;\n", in[5], in[1]);
 print("const auto tr5 = std::fma(tr4, -0.5, %s);\n", in[3]);
 print("const auto tr6 = (%24.17e) * (%s - %s);\n", sin(M_PI / 3.0), in[5], in[1]);
 print("const auto tr7 = %s + tr1;\n", in[0]);
 print("const auto tr10 = %s + tr4;\n", in[3]);
 print("%s = tr7 + tr10;\n", out[0]);
 print("%s = tr2 - tr5;\n", out[1]);
 print("%s = tr6 - tr3;\n", out[5]);
 print("%s = tr2 + tr5;\n", out[2]);
 print("%s = tr3 + tr6;\n", out[4]);
 print("%s = tr7 - tr10;\n", out[3]);
 break;
 case 7: {
 constexpr double theta = 2.0 * M_PI / 7.0;
 const double c1 = cos(theta);
 const double c2 = cos(2.0 * theta);
 const double c3 = cos(3.0 * theta);
 const double c4 = sin(theta);
 const double c5 = sin(2.0 * theta);
 const double c6 = sin(3.0 * theta);
 print("const auto tr1 = %s + %s;\n", in[1], in[6]);
 print("const auto tr2 = %s + %s;\n", in[2], in[5]);
 print("const auto tr3 = %s + %s;\n", in[3], in[4]);
 print("const auto tr4 = %s - %s;\n", in[1], in[6]);
 print("const auto tr5 = %s - %s;\n", in[2], in[5]);
 print("const auto tr6 = %s - %s;\n", in[3], in[4]);
 print("const auto tr7 = %s - 0.5 * tr3;\n", in[0]);
 print("const auto tr8 = tr1 - tr3;\n");
 print("const auto tr9 = tr2 - tr3;\n");
 print("const auto yr0 = %s + tr1 + tr2 + tr3;\n", in[0]);
 print("const auto yr1 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c1, c2);
 print("const auto yr2 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c2, c3);
 print("const auto yr3 = std::fma(tr8, %.17e, std::fma(tr9, %.17e, tr7));\n", c3, c1);
 print("const auto yr4 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c6, -c4, c5);
 print("const auto yr5 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c5, -c6, -c4);
 print("const auto yr6 = std::fma(tr4, %.17e, std::fma(tr5, %.17e, tr6 * %.17e));\n", c4, c5, c6);
 print("%s = yr0;\n", out[0]);
 print("%s = yr1;\n", out[1]);
 print("%s = -yr6;\n", out[6]);
 print("%s = yr2;\n", out[2]);
 print("%s = -yr5;\n", out[5]);
 print("%s = yr3;\n", out[3]);
 print("%s = -yr4;\n", out[4]);
 }
 break;
 case 8:
 print("const auto tr1 = %s + %s;\n", in[0], in[4]);
 print("const auto tr2 = %s - %s;\n", in[0], in[4]);
 print("const auto tr3 = %s + %s;\n", in[1], in[5]);
 print("const auto tr4 = %s - %s;\n", in[1], in[5]);
 print("const auto tr5 = %s + %s;\n", in[2], in[6]);
 print("const auto tr6 = %s - %s;\n", in[2], in[6]);
 print("const auto tr7 = %s + %s;\n", in[3], in[7]);
 print("const auto tr8 = %s - %s;\n", in[3], in[7]);
 print("const auto tr9 = tr1 + tr5;\n");
 print("const auto tr10 = tr3 + tr7;\n");
 print("const auto tr11 = (%.17e) * (tr4 - tr8);\n", 1.0 / sqrt(2));
 print("const auto tr12 = (%.17e) * (tr4 + tr8);\n", 1.0 / sqrt(2));
 print("const auto yr0 = tr9 + tr10;\n");
 print("const auto yr1 = tr2 + tr11;\n");
 print("const auto yr2 = tr1 - tr5;\n");
 print("const auto yr3 = tr2 - tr11;\n");
 print("const auto yr4 = tr9 - tr10;\n");
 print("const auto yr5 = tr12 - tr6;\n");
 print("const auto yr6 = tr3 - tr7;\n");
 print("const auto yr7 = tr12 + tr6;\n");
 print("%s = yr0;\n", out[0]);
 print("%s = yr1;\n", out[1]);
 print("%s = -yr7;\n", out[7]);
 print("%s = yr2;\n", out[2]);
 print("%s = -yr6;\n", out[6]);
 print("%s = yr3;\n", out[3]);
 print("%s = -yr5;\n", out[5]);
 print("%s = yr4;\n", out[4]);
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
 print("const auto sr1 = %s + %s;\n", in[1], in[8]);
 print("const auto sr2 = %s - %s;\n", in[1], in[8]);
 print("const auto sr3 = %s + %s;\n", in[7], in[2]);
 print("const auto sr4 = %s - %s;\n", in[7], in[2]);
 print("const auto sr5 = %s + %s;\n", in[3], in[6]);
 print("const auto sr6 = %s - %s;\n", in[3], in[6]);
 print("const auto sr7 = %s + %s;\n", in[4], in[5]);
 print("const auto sr8 = %s - %s;\n", in[4], in[5]);
 print("const auto sr9 = sr1 + sr3;\n");
 print("const auto sr10 = sr9 + sr7;\n");
 print("const auto sr11 = sr10 + sr5;\n");
 print("const auto sr12 = sr11 + %s;\n", in[0]);
 print("const auto sr13 = sr2 + sr4;\n");
 print("const auto sr14 = sr13 + sr8;\n");
 print("const auto sr15 = sr1 - sr3;\n");
 print("const auto sr16 = sr3 - sr7;\n");
 print("const auto sr17 = sr7 - sr1;\n");
 print("const auto sr18 = sr2 - sr4;\n");
 print("const auto sr19 = sr4 - sr8;\n");
 print("const auto sr20 = sr8 - sr2;\n");
 print("const auto mr0 = sr12;\n");
 print("const auto mr1 = sr10 * (%.17e);\n", ar1);
 print("const auto mi2 = sr14 * (%.17e);\n", ai2);
 print("const auto mr3 = sr5 * (%.17e);\n", ar3);
 print("const auto mi4 = sr6 * (%.17e);\n", ai4);
 print("const auto mr5 = sr15 * (%.17e);\n", ar5);
 print("const auto mr6 = sr16 * (%.17e);\n", ar6);
 print("const auto mr7 = sr17 * (%.17e);\n", ar7);
 print("const auto mi8 = sr18 * (%.17e);\n", ai8);
 print("const auto mi9 = sr19 * (%.17e);\n", ai9);
 print("const auto mi10 = sr20 * (%.17e);\n", ai10);
 print("const auto sr21 = mr1 + mr1;\n");
 print("const auto sr22 = sr21 + mr1;\n");
 print("const auto sr23 = sr22 + mr0;\n");
 print("const auto sr24 = sr23;\n");
 print("const auto si24 = mi2;\n");
 print("const auto sr25 = sr23;\n");
 print("const auto si25 = -mi2;\n");
 print("const auto sr26 = mr0 + mr3;\n");
 print("const auto sr27 = sr26 + sr21;\n");
 print("const auto sr28 = sr27 + mr5;\n");
 print("const auto sr29 = sr28 + mr6;\n");
 print("const auto sr30 = sr27 - mr6;\n");
 print("const auto sr31 = sr30 + mr7;\n");
 print("const auto sr32 = sr27 - mr5;\n");
 print("const auto sr33 = sr32 - mr7;\n");
 print("const auto si34 = mi4 + mi8;\n");
 print("const auto si35 = si34 + mi9;\n");
 print("const auto si36 = mi4 - mi9;\n");
 print("const auto si37 = si36 + mi10;\n");
 print("const auto si38 = mi4 - mi8;\n");
 print("const auto si39 = si38 - mi10;\n");
 print("const auto sr40 = sr29;\n");
 print("const auto si40 = si35;\n");
 print("const auto sr43 = sr31;\n");
 print("const auto si43 = -si37;\n");
 print("const auto sr44 = sr33;\n");
 print("const auto si44 = si39;\n");
 print("%s = mr0;\n", out[0]);
 print("%s = sr40;\n", out[1]);
 print("%s = si40;\n", out[8]);
 print("%s = sr43;\n", out[2]);
 print("%s = si43;\n", out[7]);
 print("%s = sr24;\n", out[3]);
 print("%s = si24;\n", out[6]);
 print("%s = sr44;\n", out[4]);
 print("%s = si44;\n", out[5]);

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
 print("double ");
 for (int n = 1; n < r - 2; n++) {
 print_notab("zi%i, ", n);
 }
 print_notab("zi%i;\n", r - 2);
 int N1;
 int N2;
 if (r == 10) {
 N1 = 5;
 N2 = 2;
 } else if (r == 12) {
 N1 = 3;
 N2 = 4;
 } else if (r == 14) {
 N1 = 7;
 N2 = 2;
 } else if (r == 15) {
 N1 = 5;
 N2 = 3;
 }
 std::vector<std::vector<std::string>> in1(N2, std::vector<std::string>(N1));
 std::vector<std::vector<std::string>> out1(N2, std::vector<std::string>(N1));
 std::vector<std::vector<std::string>> out2(N1 / 2 + 1, std::vector<std::string>(2 * N2));
 std::vector<std::vector<std::string>> in2(N1 / 2 + 1, std::vector<std::string>(2 * N2));
 for (int n1 = 0; n1 < N1; n1++) {
 for (int n2 = 0; n2 < N2; n2++) {
 const int nn = ((n1 * N2 + n2 * N1) % r);
 in1[n2][n1] = in[nn];
 }
 }
 for (int n2 = 0; n2 < N2; n2++) {
 out1[n2][0] = std::string("zr") + std::to_string(n2 * N1);
 for (int n1 = 1; n1 < N1 / 2 + 1; n1++) {
 const int nn = ((n1 * N2 + n2 * N1) % r);
 out1[n2][n1] = std::string("zr") + std::to_string(n2 * N1 + n1);
 if (n1 != N1 - n1) {
 out1[n2][N1 - n1] = std::string("zi") + std::to_string(n2 * N1 + n1);
 }
 }
 }
 for (int n2 = 0; n2 < N2; n2++) {
 in2[0][n2] = out1[n2][0];
 for (int n1 = 1; n1 < N1 / 2 + 1; n1++) {
 in2[n1][2 * n2] = out1[n2][n1];
 in2[n1][2 * n2 + 1] = out1[n2][N1 - n1];
 }
 }
 std::vector<std::string> negs;
 for (int n = 0; n < r; n++) {
 int n1 = n % N1;
 int n2 = n % N2;
 if (n1 == 0) {
 if (n2 > N2 / 2) {
 out2[n1][n2] = out[std::max(n, r - n)];
 negs.push_back(out[std::max(n, r - n)]);
 } else {
 out2[n1][n2] = out[std::min(n, r - n)];
 }
 } else if (n1 < N1 / 2 + 1) {
 out2[n1][2 * n2] = out[std::min(n, r - n)];
 if (n != 0) {
 out2[n1][2 * n2 + 1] = out[std::max(n, r - n)];
 if (n > r - n) {
 negs.push_back(out[n]);
 }
 }
 }

 }
 for (int n2 = 0; n2 < N2; n2++) {
 print("{\n");
 indent();
 print_real_short_fft(N1, in1[n2], out1[n2]);
 deindent();
 print("}\n");
 }
 {
 int n1 = 0;
 print("{\n");
 indent();
 print_real_short_fft(N2, in2[n1], out2[n1]);
 deindent();
 print("}\n");
 }
 for (int n1 = 1; n1 < N1 / 2 + 1; n1++) {
 print("{\n");
 indent();
 print_complex_short_fft(N2, in2[n1], out2[n1]);
 deindent();
 print("}\n");
 }
 for (int l = 0; l < negs.size(); l++) {
 print("%s = -%s;\n", negs[l].c_str(), negs[l].c_str());
 }
 }

 break;
 default:
 for (int j = 1; j <= (r - 1) / 2; j++) {
 print("const auto txp%i = %s + %s;\n", j, in[j], in[r - j]);
 }
 for (int j = 1; j <= (r - 1) / 2; j++) {
 print("const auto txm%i = %s - %s;\n", j, in[j], in[(r - j)]);
 }
 for (int i = 1; i <= (r - 1) / 2; i++) {
 print("auto ap%i = %s;\n", i, in[0]);
 print("auto bp%i = %s;\n", i, in[1]);
 for (int j = 1; j <= (r - 1) / 2; j++) {
 print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / r), i);
 if (j == 1) {
 print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / r));
 } else {
 print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / r), i);
 }
 }
 }
 print("%s = ", out[0]);
 for (int i = 0; i < r; i++) {
 print_notab(" + %s", in[i]);
 }
 print_notab(";\n");
 for (int i = 1; i <= (r - 1) / 2; i++) {
 print("%s = ap%i;\n", out[i], i);
 print("%s = -bm%i;\n", out[r - i], i);
 }
 }
 }
 */
void print_skew_short_fft(int N1, std::vector<std::string> in_, std::vector<std::string> out_) {
	print("/* begin short transform - shifted conjugate symmetric - length %i */\n", N1);
	std::vector<const char*> in(in_.size());
	std::vector<const char*> out(out_.size());
	std::vector<std::string> in1(in_.size());
	std::vector<std::string> out1(out_.size());
	for (int i = 0; i < in.size(); i++) {
		in[i] = in_[i].c_str();
		out[i] = out_[i].c_str();
	}
	if (N1 % 2 == 0 && N1 > 4) {
		std::vector<std::string> in1(N1 / 2);
		std::vector<std::string> in2(N1 / 2);
		std::vector<std::string> out1(N1 / 2);
		std::vector<std::string> out2(N1 / 2);
		print("double ar%i", 0);
		for (int n = 1; n < N1 / 2; n++) {
			print_notab(", ar%i", n);
		}
		print_notab(";\n");
		print("double ai%i", 0);
		for (int n = 1; n < N1 / 2; n++) {
			print_notab(", ai%i", n);
		}
		print_notab(";\n");
		if ((N1 / 2) % 2 == 1) {
			print("double ar%i;\n", N1 / 2);
			print("double ai%i;\n", N1 / 2);
		}
		for (int n = 0; n < N1 / 2; n++) {
			in1[n] = in_[2 * n];
			in2[n] = in_[2 * n + 1];
		}
		for (int n = 0; n < N1 / 4; n++) {
			out1[n] = std::string("ar") + std::to_string(2 * n);
			out1[N1 / 2 - n - 1] = std::string("ai") + std::to_string(2 * n);
			out2[n] = std::string("ar") + std::to_string(2 * n + 1);
			out2[N1 / 2 - n - 1] = std::string("ai") + std::to_string(2 * n + 1);
		}
		if ((N1 / 2) % 2 == 1) {
			int n = N1 / 4;
			out1[n] = std::string("ar") + std::to_string(2 * n);
			out2[n] = std::string("ar") + std::to_string(2 * n + 1);
		}
		print("{\n");
		indent();
		print_skew_short_fft(N1 / 2, in1, out1);
		deindent();
		print("}\n");
		print("{\n");
		indent();
		print_skew_short_fft(N1 / 2, in2, out2);
		deindent();
		print("}\n");
		for (int k = 0; k < N1 / 4; k++) {
			auto tw = std::polar(1.0, -2.0 * M_PI * (k + 0.5) / N1);
			int n = 2 * k + 1;
			print("tmp0 = ar%i;\n", n, n);
			print("ar%i = std::fma((%.17e), ar%i, (%.17e) * ai%i);\n", n, tw.real(), n, -tw.imag(), n);
			print("ai%i = std::fma((%.17e), ai%i, (%.17e) * tmp0);\n", n, tw.real(), n, tw.imag());
		}
		for (int n = 0; n < N1 / 4; n++) {
			print("%s = ar%i + ar%i;\n", out[n], 2 * n, 2 * n + 1);
			print("%s = ai%i + ai%i;\n", out[N1 - n - 1], 2 * n, 2 * n + 1);
		}
		int beg = N1 / 4;
		if ((N1 / 2) % 2 != 0) {
			int n = N1 / 4;
			print("%s = ar%i;\n", out[n], N1 - 2 * n - 2);
			print("%s = -ar%i;\n", out[N1 - n - 1], N1 - 2 * n - 1);
			beg++;
		}
		for (int n = beg; n < N1 / 2; n++) {
			print("%s = ar%i - ar%i;\n", out[n], N1 - 2 * n - 2, N1 - 2 * n - 1);
			print("%s = -ai%i + ai%i;\n", out[N1 - n - 1], N1 - 2 * n - 2, N1 - 2 * n - 1);
		}
		print("/* end short transform - shifted conjugate symmetric - length %i */\n", N1);
		return;
	}
	if (N1 % 2 == 1 && N1 > 4) {
		for (int n = 0; n < N1; n += 2) {
			print("const auto a%i = %s;\n", n / 2, in[n]);
			in1[n / 2] = std::string("a") + std::to_string(n / 2);
		}
		for (int n = 1; n < N1; n += 2) {
			print("const auto a%i = -%s;\n", (n + N1) / 2, in[n]);
			in1[(n + N1) / 2] = std::string("a") + std::to_string((n + N1) / 2);
		}
		for (int n = 0; n < (N1 + 1) / 4; n++) {
			out1[2 * n + 1] = out[n];
			out1[N1 - 2 * n - 1] = out[N1 - n - 1];
		}
		for (int n = 0; n < (N1) / 4; n++) {
			out1[2 * (N1 / 4 - n)] = out[n + (N1 + 1) / 4];
			out1[N1 - 2 * (N1 / 4 - n)] = out[N1 - n - 1 - (N1 + 1) / 4];
		}
		for (int n = 0; n < (N1 + 1) / 4; n++) {
		}
		out1[0] = out[N1 / 2];
		print("{\n");
		indent();
		print_real_short_fft(N1, in1, out1);
		deindent();
		print("}\n");
		for (int n = 0; n < (N1) / 4; n++) {
			print("%s = -%s;\n", out[N1 - n - 1 - (N1 + 1) / 4], out[N1 - n - 1 - (N1 + 1) / 4]);
		}
		print("/* end short transform - shifted conjugate symmetric - length %i */\n", N1);
		return;
	}

	switch (N1) {
	case 1:
		return;
	case 2:
		for (int k1 = 0; k1 < N1; k1++) {
			print("const auto x%i = %s;\n", k1, in[k1]);
		}
		print("%s = x0;\n", out[0]);
		print("%s = -x1;\n", out[1]);
		break;
	case 3:
		for (int k1 = 0; k1 < N1; k1++) {
			print("const auto x%i = %s;\n", k1, in[k1]);
		}
		print("const auto t1 = x1 - x2;\n");
		print("%s = x0 + 0.5 * t1;\n", out[0]);
		print("%s = x0 - t1;\n", out[2]);
		print("%s = (%.17e) * (x1 + x2);\n", out[1]);
		break;
	case 4:
		for (int k1 = 0; k1 < N1; k1++) {
			print("const auto x%i = %s;\n", k1, in[k1]);
		}
		print("const auto t1 = (x1 - x3) * (%.17e);\n", 1.0 / sqrt(2));
		print("const auto t2 = (x1 + x3) * (%.17e);\n", 1.0 / sqrt(2));
		print("%s = x0 + t1;\n", out[0]);
		print("%s = -x2 - t2;\n", out[3]);
		print("%s = x0 - t1;\n", out[1]);
		print("%s = x2 - t2;\n", out[2]);
		break;
	case 5:
		for (int k1 = 0; k1 < N1; k1++) {
			print("const auto x%i = %s;\n", k1, in[k1]);
		}
		print("const auto t1 = x1 - x4;\n");
		print("const auto t2 = x1 + x4;\n");
		print("const auto t3 = x2 - x3;\n");
		print("const auto t4 = x2 + x3;\n");
		print("const auto t5 = t1 - t3;\n");
		print("const auto t6 = std::fma(0.25, t5, x0);\n");
		print("const auto t7 = (%.17e) * (t1 + t3);\n", sqrt(5) * 0.25);
		print("%s = t6 + t7;\n", out[0]);
		print("%s = std::fma((%.17e), t2, (%.17e) * t4);\n", out[4], -sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
		print("%s = t6 - t7;\n", out[1]);
		print("%s = std::fma((%.17e), t2, (%.17e) * t4);\n", out[3], -sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
		print("%s = x0 - t5;\n", out[3]);
		break;
	case 6:
		for (int k1 = 0; k1 < N1; k1++) {
			print("const auto x%i = %s;\n", k1, in[k1]);
		}
		print("const auto t1 = (%.17e) * (x5 - x1);\n", sin(M_PI / 3.0));
		print("const auto t2 = (%.17e) * (x2 + x4);\n", sin(M_PI / 3.0));
		print("const auto t3 = x2 - x4;\n");
		print("const auto t4 = x1 + x5;\n");
		print("const auto t5 = std::fma(0.5, t3, x0);\n");
		print("const auto t6 = -std::fma(0.5, t4, x3);\n");
		print("%s = t5 - t1;\n", out[0]);
		print("%s = t6 - t2;\n", out[5]);
		print("%s = x0 - t3;\n", out[1]);
		print("%s = x3 - t4;\n", out[4]);
		print("%s = t5 + t1;\n", out[2]);
		print("%s = t6 + t2;\n", out[3]);
		break;
	default:
		abort();

	}
	print("/* end short transform - shifted conjugate symmetric - length %i */\n", N1);
}

