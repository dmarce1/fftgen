#include "fftgen.hpp"

#include <set>

bool use_fma = true;

void print_vars(std::vector<std::string> vars) {
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
}

int print_short_fft(int N, std::vector<std::string> in_, std::vector<std::string> out_, bool real = false) {
	int ocnt = 0;
	std::vector<const char*> in(in_.size());
	std::vector<const char*> out(out_.size());
	for (int i = 0; i < in.size(); i++) {
		in[i] = in_[i].c_str();
		out[i] = out_[i].c_str();
	}
	if (N == 2) {
		print("/* begin short transform - %s - length %i */\n", real ? "real" : "complex", N);
		if (!real) {
			print("const auto tr = %s;\n", in[2]);
			print("const auto ti = %s;\n", in[3]);
			print("%s = %s - tr;\n", out[2], in[0]);
			print("%s = %s - ti;\n", out[3], in[1]);
			print("%s = %s + tr;\n", out[0], in[0]);
			print("%s = %s + ti;\n", out[1], in[1]);
			print("/* end short transform */\n");
			return 4;
		} else {
			print("const auto t = %s;\n", in[1]);
			print("%s = %s - t;\n", out[1], in[0]);
			print("%s = %s + t;\n", out[0], in[0]);
			print("/* end short transform */\n");
			return 2;
		}
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
		print("/* begin short transform - %s - length %i */\n", real ? "real" : "complex", N);
		if (pfac[0].second > 1 /*&& !(real && N % 2 == 0)*/) {
			if (!real && N % 4 == 0 && N >= 4) {
				std::vector<std::string> in0;
				std::vector<std::string> in1;
				std::vector<std::string> in3;
				std::vector<std::string> out0;
				std::vector<std::string> out1;
				std::vector<std::string> out3;
				std::vector<std::string> vars;
				for (int n = 0; n < N; n++) {
					vars.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(n));
					vars.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(n));
				}
				print_vars(vars);
				for (int n = 0; n < N; n += 2) {
					in0.push_back(in[2 * n]);
					in0.push_back(in[2 * n + 1]);
					out0.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(n));
					out0.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(n));
				}
				for (int n = 1; n < N; n += 4) {
					in1.push_back(in[2 * n]);
					in1.push_back(in[2 * n + 1]);
					out1.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(n));
					out1.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(n));
				}
				for (int n = 3; n < N; n += 4) {
					in3.push_back(in[2 * n]);
					in3.push_back(in[2 * n + 1]);
					out3.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(n));
					out3.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(n));
				}
				if (N > 4) {
					print("{\n");
					indent();
					ocnt += print_short_fft(N / 4, in1, out1);
					deindent();
					print("}\n");
				} else {
					print("%s = %s;\n", out1[0].c_str(), in1[0].c_str());
					print("%s = %s;\n", out1[1].c_str(), in1[1].c_str());
				}
				print("{\n");
				indent();
				ocnt += print_short_fft(N / 2, in0, out0);
				deindent();
				print("}\n");
				if (N > 4) {
					print("{\n");
					indent();
					ocnt += print_short_fft(N / 4, in3, out3);
					deindent();
					print("}\n");
				} else {
					print("%s = %s;\n", out3[0].c_str(), in3[0].c_str());
					print("%s = %s;\n", out3[1].c_str(), in3[1].c_str());
				}
				for (int n = 1; n < N / 4; n++) {
					for (int k = 1; k < 4; k += 2) {
						const auto tw = std::polar(1.0, -2.0 * M_PI * n * k / N);
						int q = 4 * n + k;
						print("tmp0 = y%ir%i;\n", N, q);
						print("y%ir%i = std::fma(y%ir%i, %.17e, y%ii%i * %.17e);\n", N, q, N, q, tw.real(), N, q, -tw.imag());
						print("y%ii%i = std::fma(tmp0, %.17e, y%ii%i * %.17e);\n", N, q, tw.imag(), N, q, tw.real());
						ocnt += (6 - 2 * use_fma);
					}
				}
				for (int k = 0; k < N / 4; k++) {
					print("{\n");
					indent();
					print("const auto zrs = (y%ir%i + y%ir%i);\n", N, 4 * k + 1, N, 4 * k + 3);
					print("const auto zis = (y%ii%i + y%ii%i);\n", N, 4 * k + 1, N, 4 * k + 3);
					print("const auto zrd = (y%ir%i - y%ir%i);\n", N, 4 * k + 1, N, 4 * k + 3);
					print("const auto zid = (y%ii%i - y%ii%i);\n", N, 4 * k + 1, N, 4 * k + 3);
					print("%s = y%ir%i + zrs;\n", out[2 * k + 0 * N / 2 + 0], N, 2 * k + 0 * N / 2);
					print("%s = y%ii%i + zis;\n", out[2 * k + 0 * N / 2 + 1], N, 2 * k + 0 * N / 2);
					print("%s = y%ir%i + zid;\n", out[2 * k + 1 * N / 2 + 0], N, 2 * k + 1 * N / 2);
					print("%s = y%ii%i - zrd;\n", out[2 * k + 1 * N / 2 + 1], N, 2 * k + 1 * N / 2);
					print("%s = y%ir%i - zrs;\n", out[2 * k + 2 * N / 2 + 0], N, 2 * k + 0 * N / 2);
					print("%s = y%ii%i - zis;\n", out[2 * k + 2 * N / 2 + 1], N, 2 * k + 0 * N / 2);
					print("%s = y%ir%i - zid;\n", out[2 * k + 3 * N / 2 + 0], N, 2 * k + 1 * N / 2);
					print("%s = y%ii%i + zrd;\n", out[2 * k + 3 * N / 2 + 1], N, 2 * k + 1 * N / 2);
					deindent();
					print("}\n");
					ocnt += 12;
				}
			} else if (!real) {
				std::vector<std::string> vars;
				int R = pfac[0].first;
				std::vector<std::vector<std::string>> in1(R);
				std::vector<std::vector<std::string>> out1(R);
				for (int n = 0; n < N; n++) {
					vars.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(n));
					vars.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(n));
				}
				print_vars(vars);
				for (int n = 0; n < N / R; n++) {
					for (int r = 0; r < R; r++) {
						in1[r].push_back(in[2 * (R * n + r)]);
						in1[r].push_back(in[2 * (R * n + r) + 1]);
						out1[r].push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * n + r));
						out1[r].push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(R * n + r));
					}
				}
				if (N > R) {
					for (int r = 0; r < R; r++) {
						print("{\n");
						indent();
						ocnt += print_short_fft(N / R, in1[r], out1[r]);
						deindent();
						print("}\n");
					}
				}
				for (int n = 1; n < N / R; n++) {
					for (int k = 1; k < R; k++) {
						const auto tw = std::polar(1.0, -2.0 * M_PI * n * k / N);
						int q = R * n + k;
						if (n * k == 0) {

						} else if (n * k == N / 2 && N % 2 == 0) {
							print("y%ir%i = -y%ir%i;\n", N, q, N, q);
							print("y%ii%i = -y%ii%i;\n", N, q, N, q);
							ocnt += 2;
						} else if (n * k == N / 4 && N % 4 == 0) {
							print("std::swap(y%ir%i, y%ii%i);\n", N, q, N, q);
							print("y%ii%i = -y%ii%i;\n", N, q, N, q);
							ocnt++;
						} else if (n * k == 3 * N / 4 && N % 4 == 0) {
							print("std::swap(y%ir%i, y%ii%i);\n", N, q, N, q);
							print("y%ir%i = -y%ir%i;\n", N, q, N, q);
							ocnt++;
						} else {
							print("tmp0 = y%ir%i;\n", N, q);
							print("y%ir%i = std::fma(y%ir%i, %.17e, y%ii%i * %.17e);\n", N, q, N, q, tw.real(), N, q, -tw.imag());
							print("y%ii%i = std::fma(tmp0, %.17e, y%ii%i * %.17e);\n", N, q, tw.imag(), N, q, tw.real());
							ocnt += (6 - 2 * use_fma);
						}
					}
				}
				for (int k = 0; k < N / R; k++) {
					print("{\n");
					indent();
					std::vector<std::string> in2, out2;
					for (int r = 0; r < R; r++) {
						in2.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * k + r));
						in2.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(R * k + r));
						out2.push_back(out[2 * (N / R * r + k)]);
						out2.push_back(out[2 * (N / R * r + k) + 1]);
					}
					ocnt += print_short_fft(R, in2, out2);
					deindent();
					print("}\n");
				}
			} else {
				if (N % 4 == 0) {
					std::vector<std::string> vars;
					int R = 4;
					std::vector<std::vector<std::string>> in1(R);
					std::vector<std::vector<std::string>> out1(R);
					for (int n = 0; n < N; n++) {
						vars.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(n));
					}
					for (int n = 0; n < N; n++) {
						vars.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(n));
					}
					print_vars(vars);
					for (int n = 0; n < N / R; n++) {
						for (int r = 0; r < R; r++) {
							in1[r].push_back(in[R * n + r]);
							if (n < N / R / 2 + 1) {
								out1[r].push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * n + r));
							} else {
								out1[r].push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(R * (N / R - n) + r));
							}
						}
					}
					if (N > R) {
						for (int r = 0; r < R; r++) {
							print("{\n");
							indent();
							ocnt += print_real_short_fft(N / R, in1[r], out1[r]);
							deindent();
							print("}\n");
						}
					}
					for (int k = 1; k < R; k++) {
						for (int n = 1; n <= ((N / R) - 1) / 2; n++) {
							const auto tw = std::polar(1.0, -2.0 * M_PI * n * k / N);
							int q = R * n + k;
							if (n * k == 0) {

							} else if (n * k == N / 2 && N % 2 == 0) {
								print("y%ir%i = -y%ir%i;\n", N, q, N, q);
								print("y%ii%i = -y%ii%i;\n", N, q, N, q);
								ocnt += 2;
							} else if (n * k == N / 4 && N % 4 == 0) {
								print("std::swap(y%ir%i, y%ii%i);\n", N, q, N, q);
								print("y%ii%i = -y%ii%i;\n", N, q, N, q);
								ocnt++;
							} else if (n * k == 3 * N / 4 && N % 4 == 0) {
								print("std::swap(y%ir%i, y%ii%i);\n", N, q, N, q);
								print("y%ir%i = -y%ir%i;\n", N, q, N, q);
								ocnt++;
							} else {
								print("tmp0 = y%ir%i;\n", N, q);
								print("y%ir%i = std::fma(y%ir%i, %.17e, y%ii%i * %.17e);\n", N, q, N, q, tw.real(), N, q, -tw.imag());
								print("y%ii%i = std::fma(tmp0, %.17e, y%ii%i * %.17e);\n", N, q, tw.imag(), N, q, tw.real());
								ocnt += (6 - 2 * use_fma);
							}
						}
					}
					{
						int k = 0;
						print("{\n");
						indent();
						std::vector<std::string> in2, out2;
						for (int r = 0; r < R; r++) {
							in2.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * k + r));
							out2.push_back(out[(N / R * r + k)]);
						}
						ocnt += print_real_short_fft(R, in2, out2);
						deindent();
						print("}\n");
					}
					std::vector<std::string> negs;
					for (int k = 1; k <= ((N / R) - 1) / 2; k++) {
						print("{\n");
						indent();
						std::vector<std::string> in2, out2;
						for (int r = 0; r < R; r++) {
							in2.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * k + r));
							in2.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(R * k + r));
							int iii = (N / R * r + k);
							if (iii < N / 2 + 1) {
								out2.push_back(out[iii]);
								out2.push_back(out[N - iii]);
							} else {
								out2.push_back(out[N - iii]);
								out2.push_back(out[iii]);
								negs.push_back(out[iii]);
							}

						}
						ocnt += print_short_fft(R, in2, out2);
						deindent();
						print("}\n");
					}
					if ((N / R) % 2 == 0) {
						int k = N / R / 2;
						print("{\n");
						indent();
						std::vector<std::string> in2, out2;
						for (int r = 0; r < R; r++) {
							in2.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * k + r));
							out2.push_back(out[(N / R * r + k)]);
						}
						ocnt += print_skew_short_fft(R, in2, out2);
						deindent();
						print("}\n");
					}
					for (int i = 0; i < negs.size(); i++) {
						print("%s = -%s;\n", negs[i].c_str(), negs[i].c_str());
					}
				} else {
					std::vector<std::string> vars;
					int R = pfac[0].first;
					std::vector<std::vector<std::string>> in1(R);
					std::vector<std::vector<std::string>> out1(R);
					for (int n = 0; n < N; n++) {
						vars.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(n));
					}
					for (int n = 0; n < N; n++) {
						vars.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(n));
					}
					print_vars(vars);
					for (int n = 0; n < N / R; n++) {
						for (int r = 0; r < R; r++) {
							in1[r].push_back(in[R * n + r]);
							if (n < N / R / 2 + 1) {
								out1[r].push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * n + r));
							} else {
								out1[r].push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(R * (N / R - n) + r));
							}
						}
					}
					if (N > R) {
						for (int r = 0; r < R; r++) {
							print("{\n");
							indent();
							ocnt += print_real_short_fft(N / R, in1[r], out1[r]);
							deindent();
							print("}\n");
						}
					}
					for (int k = 1; k < R; k++) {
						for (int n = 1; n <= ((N / R) - 1) / 2; n++) {
							const auto tw = std::polar(1.0, -2.0 * M_PI * n * k / N);
							int q = R * n + k;
							if (n * k == 0) {

							} else if (n * k == N / 2 && N % 2 == 0) {
								print("y%ir%i = -y%ir%i;\n", N, q, N, q);
								print("y%ii%i = -y%ii%i;\n", N, q, N, q);
								ocnt += 2;
							} else if (n * k == N / 4 && N % 4 == 0) {
								print("std::swap(y%ir%i, y%ii%i);\n", N, q, N, q);
								print("y%ii%i = -y%ii%i;\n", N, q, N, q);
								ocnt++;
							} else if (n * k == 3 * N / 4 && N % 4 == 0) {
								print("std::swap(y%ir%i, y%ii%i);\n", N, q, N, q);
								print("y%ir%i = -y%ir%i;\n", N, q, N, q);
								ocnt++;
							} else {
								print("tmp0 = y%ir%i;\n", N, q);
								print("y%ir%i = std::fma(y%ir%i, %.17e, y%ii%i * %.17e);\n", N, q, N, q, tw.real(), N, q, -tw.imag());
								print("y%ii%i = std::fma(tmp0, %.17e, y%ii%i * %.17e);\n", N, q, tw.imag(), N, q, tw.real());
								ocnt += (6 - 2 * use_fma);
							}
						}
					}
					{
						int k = 0;
						print("{\n");
						indent();
						std::vector<std::string> in2, out2;
						for (int r = 0; r < R; r++) {
							in2.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * k + r));
							out2.push_back(out[(N / R * r + k)]);
						}
						ocnt += print_real_short_fft(R, in2, out2);
						deindent();
						print("}\n");
					}
					std::vector<std::string> negs;
					for (int k = 1; k <= ((N / R) - 1) / 2; k++) {
						print("{\n");
						indent();
						std::vector<std::string> in2, out2;
						for (int r = 0; r < R; r++) {
							in2.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * k + r));
							in2.push_back(std::string("y") + std::to_string(N) + "i" + std::to_string(R * k + r));
							int iii = (N / R * r + k);
							if (iii < N / 2 + 1) {
								out2.push_back(out[iii]);
								out2.push_back(out[N - iii]);
							} else {
								out2.push_back(out[N - iii]);
								out2.push_back(out[iii]);
								negs.push_back(out[iii]);
							}

						}
						ocnt += print_short_fft(R, in2, out2);
						deindent();
						print("}\n");
					}
					if ((N / R) % 2 == 0) {
						int k = N / R / 2;
						print("{\n");
						indent();
						std::vector<std::string> in2, out2;
						for (int r = 0; r < R; r++) {
							in2.push_back(std::string("y") + std::to_string(N) + "r" + std::to_string(R * k + r));
							out2.push_back(out[(N / R * r + k)]);
						}
						ocnt += print_skew_short_fft(R, in2, out2);
						deindent();
						print("}\n");
					}
					for (int i = 0; i < negs.size(); i++) {
						print("%s = -%s;\n", negs[i].c_str(), negs[i].c_str());
					}

				}
			}
		} else {

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
			print_vars(vars);
			if (!real) {
				for (int j = 1; j <= (N - 1) / 2; j++) {
					print("const auto tp%i = %s + %s;\n", j, in[ri[j]], in[ri[N - j]]);
					ocnt++;
				}
				for (int j = 1; j <= (N - 1) / 2; j++) {
					print("const auto tm%i = %s - %s;\n", j, in[ri[j]], in[ri[N - j]]);
					ocnt++;
				}
				for (int j = 1; j <= (N - 1) / 2; j++) {
					print("const auto up%i = %s + %s;\n", j, in[ii[j]], in[ii[N - j]]);
					ocnt++;
				}
				for (int j = 1; j <= (N - 1) / 2; j++) {
					print("const auto um%i = %s - %s;\n", j, in[ii[j]], in[ii[N - j]]);
					ocnt++;
				}
			} else {
				for (int j = 1; j <= (N - 1) / 2; j++) {
					print("const auto tp%i = %s + %s;\n", j, in[ri[j]], in[ii[N - j]]);
					ocnt++;
				}
				for (int j = 1; j <= (N - 1) / 2; j++) {
					ocnt++;
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
						ocnt += real ? 1 : 2;
						real ? void() : print("bp%i += up%i;\n", i, j);
					} else if (k == N / 2 && N % 2 == 0) {
						ocnt += real ? 1 : 2;
						print("ap%i -= tp%i;\n", i, j);
						real ? void() : print("bp%i -= up%i;\n", i, j);
					} else if (k == N / 4 && N % 4 == 0) {
					} else if (k == 3 * N / 4 && N % 4 == 0) {
					} else {
						ocnt += real ? (2 - use_fma) : (4 - 2 * use_fma);
						print("ap%i = std::fma(tp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / N), i);
						real ? void() : print("bp%i = std::fma(up%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / N), i);
					}
				}
				if (N % 2 == 0) {
					print("ap%i %c= %s;\n", i, i % 2 == 0 ? '+' : '-', in[ri[N / 2]]);
					ocnt += real ? 1 : 2;
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
							ocnt += real ? 1 : 2;
							init = true;
							init_ele[j] = true;
							break;
						} else {
							real ? void() : print("am%i = um%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N));
							print("bm%i = tm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N));
							init = true;
							ocnt += real ? 1 : 2;
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
							ocnt += real ? 1 : 2;
						} else if (k == 3 * N / 4 && N % 4 == 0) {
							real ? void() : print("am%i -= um%i;\n", i, j);
							ocnt += real ? 1 : 2;
							print("bm%i -= tm%i;\n", i, j);
						} else {
							real ? void() : print("am%i = std::fma(um%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / N), i);
							ocnt += real ? (2 - use_fma) : (4 - 2 * use_fma);
							print("bm%i = std::fma(tm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / N), i);
						}
					}
				}
			}
			if (N % 2 == 0) {
				if (N % 4 == 0) {
					print("pr0 = sr0 = %s + %s;\n", in[ri[0]], in[ri[N / 2]]);
					ocnt++;
					if (!real) {
						print("pi0 = si0 = %s + %s;\n", in[ii[0]], in[ii[N / 2]]);
						ocnt++;
					}
				} else {
					print("sr0 = %s + %s;\n", in[ri[0]], in[ri[N / 2]]);
					ocnt++;
					print("pr0 = %s - %s;\n", in[ri[0]], in[ri[N / 2]]);
					ocnt++;
					if (!real) {
						print("si0 = %s + %s;\n", in[ii[0]], in[ii[N / 2]]);
						ocnt++;
						print("pi0 = %s - %s;\n", in[ii[0]], in[ii[N / 2]]);
						ocnt++;
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
					ocnt++;
				}
				for (int i = 3; i <= (N - 1) / 2; i += 2) {
					print("or0 += tp%i;\n", i);
					ocnt++;
				}
				print("sr0 += er0 + or0;\n");
				ocnt++;
				ocnt++;
				print("pr0 += er0 - or0;\n");
				ocnt++;
				ocnt++;
				if (!real) {
					for (int i = 4; i <= (N - 1) / 2; i += 2) {
						print("ei0 += up%i;\n", i);
						ocnt++;
					}
					for (int i = 3; i <= (N - 1) / 2; i += 2) {
						print("oi0 += up%i;\n", i);
						ocnt++;
					}
					print("si0 += ei0 + oi0;\n");
					print("pi0 += ei0 - oi0;\n");
					ocnt++;
					ocnt++;
					ocnt++;
					ocnt++;
				}
			} else {
				for (int i = 1; i <= (N - 1) / 2; i++) {
					print("sr0 += tp%i;\n", i);
					ocnt++;
				}
				if (!real) {
					for (int i = 1; i <= (N - 1) / 2; i++) {
						print("si0 += up%i;\n", i);
						ocnt++;
					}
				}
				if (N % 2 == 0) {
					for (int i = 1; i <= (N - 1) / 2; i++) {
						print("pr0 %c= tp%i;\n", i % 2 == 0 ? '+' : '-', i);
						ocnt++;
					}
					if (!real) {
						for (int i = 1; i <= (N - 1) / 2; i++) {
							print("pi0 %c= up%i;\n", i % 2 == 0 ? '+' : '-', i);
							ocnt++;
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
					ocnt += 4;
				}
			} else {
				print("%s = sr0;\n", out[ri[0]]);
				if (N % 2 == 0) {
					print("%s = pr0;\n", out[ri[N / 2]]);
				}
				for (int i = 1; i <= (N - 1) / 2; i++) {
					print("%s = ap%i;\n", out[ri[i]], i, i);
					print("%s = -bm%i;\n", out[ii[i]], i, i);
					ocnt++;
				}
			}
			//	}
		}
		print("/* end short transform */\n");
	} else {
		if (real) {
			int N2 = pow(pfac[0].first, pfac[0].second);
			int N1 = N / N2;
			print("/* begin compound short transform - %s - length %i = %i x %i */\n", real ? "real" : "complex", N, N1, N2);
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
				int sgn = 1;
				if (n1 < N1 / 2 + 1) {
					if (n1 == 0) {
						out2[n1][std::min(n2, N2 - n2)] = out[std::min(n, N - n)];
						if (n != 0) {
							out2[n1][std::max(n2, N2 - n2)] = out[std::max(n, N - n)];
						}
						if (n2 > N2 - n2 && n <= N - n) {
							negs.push_back(out[std::max(n, N - n)]);
						}
					} else {
						out2[n1][2 * n2] = out[std::min(n, N - n)];
						if (n != 0) {
							out2[n1][2 * n2 + 1] = out[std::max(n, N - n)];
						}
						if (n > N - n) {
							negs.push_back(out[std::max(n, N - n)]);
						}
					}
				}

			}
			for (int n2 = 0; n2 < N2; n2++) {
				print("{\n");
				indent();
				ocnt += print_real_short_fft(N1, in1[n2], out1[n2]);
				deindent();
				print("}\n");
			}
			{
				int n1 = 0;
				print("{\n");
				indent();
				ocnt += print_real_short_fft(N2, in2[n1], out2[n1]);
				deindent();
				print("}\n");
			}
			for (int n1 = 1; n1 < N1 / 2 + 1; n1++) {
				print("{\n");
				indent();
				ocnt += print_complex_short_fft(N2, in2[n1], out2[n1]);
				deindent();
				print("}\n");
			}
			for (int l = 0; l < negs.size(); l++) {
				print("%s = -%s;\n", negs[l].c_str(), negs[l].c_str());
				ocnt++;
			}
			print("/* end short transform */\n");
		} else {
			int N1 = pow(pfac[0].first, pfac[0].second);
			int N2 = N / N1;
			print("/* begin compound short transform - %s - length %i = %i x %i */\n", real ? "real" : "complex", N, N1, N2);

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
				ocnt += print_complex_short_fft(N1, in1[n2], out1[n2]);
				deindent();
				print("}\n");
			}
			for (int n1 = 0; n1 < N1; n1++) {
				print("{\n");
				indent();
				ocnt += print_complex_short_fft(N2, in2[n1], out2[n1]);
				deindent();
				print("}\n");
			}
			print("/* end short transform */\n");
		}
	}
	static std::set<int> rdone;
	static std::set<int> cdone;
	if (!real) {
		if (cdone.find(N) == cdone.end()) {
			cdone.insert(N);
			fprintf( stderr, "COMPLEX OP COUNT - %i = %i\n", N, ocnt);
		}
	} else {
		if (rdone.find(N) == rdone.end()) {
			rdone.insert(N);
			fprintf( stderr, "  REAL  OP COUNT - %i = %i\n", N, ocnt);
		}
	}
	return ocnt;
}

int print_complex_short_fft(int N, std::vector<std::string> in_, std::vector<std::string> out_) {
	return print_short_fft(N, in_, out_);
}

int print_real_short_fft(int N, std::vector<std::string> in_, std::vector<std::string> out_) {
	return print_short_fft(N, in_, out_, true);
}

int print_skew_short_fft(int N1, std::vector<std::string> in_, std::vector<std::string> out_) {
	int ocnt = 0;
	print("/* begin short transform - shifted conjugate symmetric - length %i */\n", N1);
	std::vector<const char*> in(in_.size());
	std::vector<const char*> out(out_.size());
	std::vector<std::string> in1(in_.size());
	std::vector<std::string> out1(out_.size());
	for (int i = 0; i < in.size(); i++) {
		in[i] = in_[i].c_str();
		out[i] = out_[i].c_str();
	}
	if (N1 == 2) {
		for (int k1 = 0; k1 < N1; k1++) {
			print("const auto x%i = %s;\n", k1, in[k1]);
		}
		print("%s = x0;\n", out[0]);
		print("%s = -x1;\n", out[1]);
		ocnt++;
	} else if (N1 % 2 == 0) {
		int N = N1;
		std::vector<std::string> in1(N1 / 2);
		std::vector<std::string> in2(N1 / 2);
		std::vector<std::string> out1(N1 / 2);
		std::vector<std::string> out2(N1 / 2);
		print("double a%ir%i", N, 0);
		for (int n = 1; n < N1 / 2; n++) {
			print_notab(", a%ir%i", N, n);
		}
		print_notab(";\n");
		print("double a%ii%i", N, 0);
		for (int n = 1; n < N1 / 2; n++) {
			print_notab(", a%ii%i", N, n);
		}
		print_notab(";\n");
		if ((N1 / 2) % 2 == 1) {
			print("double a%ir%i;\n", N, N1 / 2);
			print("double a%ii%i;\n", N, N1 / 2);
		}
		for (int n = 0; n < N1 / 2; n++) {
			in1[n] = in_[2 * n];
			in2[n] = in_[2 * n + 1];
		}
		for (int n = 0; n < N1 / 4; n++) {
			out1[n] = std::string("a") + std::to_string(N) + "r" + std::to_string(2 * n);
			out1[N1 / 2 - n - 1] = std::string("a") + std::to_string(N) + "i" + std::to_string(2 * n);
			out2[n] = std::string("a") + std::to_string(N) + "r" + std::to_string(2 * n + 1);
			out2[N1 / 2 - n - 1] = std::string("a") + std::to_string(N) + "i" + std::to_string(2 * n + 1);
		}
		if ((N1 / 2) % 2 == 1) {
			int n = N1 / 4;
			out1[n] = std::string("a") + std::to_string(N) + "r" + std::to_string(2 * n);
			out2[n] = std::string("a") + std::to_string(N) + "r" + std::to_string(2 * n + 1);
		}
		print("{\n");
		indent();
		ocnt += print_skew_short_fft(N1 / 2, in1, out1);
		deindent();
		print("}\n");
		print("{\n");
		indent();
		ocnt += print_skew_short_fft(N1 / 2, in2, out2);
		deindent();
		print("}\n");
		for (int k = 0; k < N1 / 4; k++) {
			auto tw = std::polar(1.0, -2.0 * M_PI * (k + 0.5) / N1);
			int n = 2 * k + 1;
			print("tmp0 = a%ir%i;\n", N, n);
			print("a%ir%i = std::fma((%.17e), a%ir%i, (%.17e) * a%ii%i);\n", N, n, tw.real(), N, n, -tw.imag(), N, n);
			print("a%ii%i = std::fma((%.17e), a%ii%i, (%.17e) * tmp0);\n", N, n, tw.real(), N, n, tw.imag());
			ocnt += 6 - 2 * use_fma;
		}
		for (int n = 0; n < N1 / 4; n++) {
			print("%s = a%ir%i + a%ir%i;\n", out[n], N, 2 * n, N, 2 * n + 1);
			print("%s = a%ii%i + a%ii%i;\n", out[N1 - n - 1], N, 2 * n, N, 2 * n + 1);
			ocnt += 2;
		}
		int beg = N1 / 4;
		if ((N1 / 2) % 2 != 0) {
			int n = N1 / 4;
			print("%s = a%ir%i;\n", out[n], N, N1 - 2 * n - 2);
			print("%s = -a%ir%i;\n", out[N1 - n - 1], N, N1 - 2 * n - 1);
			ocnt++;
			beg++;
		}
		for (int n = beg; n < N1 / 2; n++) {
			print("%s = a%ir%i - a%ir%i;\n", out[n], N, N1 - 2 * n - 2, N, N1 - 2 * n - 1);
			print("%s = a%ii%i - a%ii%i;\n", out[N1 - n - 1], N, N1 - 2 * n - 1, N, N1 - 2 * n - 2);
			ocnt += 2;
		}
	} else {
		for (int n = 0; n < N1; n += 2) {
			print("const auto a%i = %s;\n", n / 2, in[n]);
			in1[n / 2] = std::string("a") + std::to_string(n / 2);
		}
		for (int n = 1; n < N1; n += 2) {
			print("const auto a%i = -%s;\n", (n + N1) / 2, in[n]);
			in1[(n + N1) / 2] = std::string("a") + std::to_string((n + N1) / 2);
			ocnt++;
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
		print("/* end short transform */\n");
	}
	return ocnt;
}

