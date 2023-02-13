#include "fftgen.hpp"

void fft_real(int N) {
	int N1;
	N1 = 1;
	while (!(N % N1 == 0 && is_prime(N1))) {
		N1++;
	}
/*	if (N1 < 6 && N % 6 == 0) {
	//	N1 = 6;
	}
	if (N1 < 4 && N % 4 == 0) {
		N1 = 4;
	}*/

	int N2 = N / N1;
	print("{\n");
	indent();

	print("std::array<std::array<double, %i>, %i> z0;\n", 2 * N2, N1 / 2);
	if (N1 % 2 != 0) {
		print("std::array<double, %i> z1;\n", N2);
	}
	for (int n2 = 0; n2 < N2; n2++) {
		for (int n1 = 0; n1 < N1; n1 += 2) {
			if (N1 - n1 == 1) {
				break;
			}
			print("z0[%i][%i] = x[%i];\n", n1 / 2, 2 * n2, N1 * n2 + n1);
			print("z0[%i][%i] = x[%i];\n", n1 / 2, 2 * n2 + 1, N1 * n2 + n1 + 1);
		}
		if (N1 % 2 != 0) {
			const int n1 = N1 - 1;
			print("z1[%i] = x[%i];\n", n2, N1 * n2 + n1);
		}
	}
	if (N2 > 1) {
		for (int n1 = 0; n1 < N1; n1 += 2) {
			if (N1 - n1 == 1) {
				break;
			}
			print("fft_%i(z0[%i].data());\n", N2, n1 / 2);
		}
		if (N1 % 2 != 0) {
			print("fft_real_base_%i(z1.data());\n", N2);
		}
	}
	for (int k2 = 0; k2 < N2 / 2 + 1; k2++) {
		for (int n1 = 0; n1 < N1; n1 += 2) {
			const int nkr0 = N2 * n1 + k2;
			const int nkr1 = N2 * (n1 + 1) + k2;
			const int nki0 = N2 * n1 + (N2 - k2);
			const int nki1 = N2 * (n1 + 1) + (N2 - k2);
			if (N1 - n1 == 1) {
				break;
			}
			if (k2 == 0 || (k2 == N2 / 2 && N2 % 2 == 0)) {
				print("x[%i] = z0[%i][%i];\n", nkr0, n1 / 2, 2 * k2);
				print("x[%i] = z0[%i][%i];\n", nkr1, n1 / 2, 2 * k2 + 1);
			} else {
				print("x[%i] = 0.5 * (z0[%i][%i] + z0[%i][%i]);\n", nkr0, n1 / 2, 2 * k2, n1 / 2, 2 * (N2 - k2));
				print("x[%i] = 0.5 * (z0[%i][%i] - z0[%i][%i]);\n", nki0, n1 / 2, 2 * k2 + 1, n1 / 2, 2 * (N2 - k2) + 1);
				print("x[%i] = 0.5 * (z0[%i][%i] + z0[%i][%i]);\n", nkr1, n1 / 2, 2 * k2 + 1, n1 / 2, 2 * (N2 - k2) + 1);
				print("x[%i] = 0.5 * (z0[%i][%i] - z0[%i][%i]);\n", nki1, n1 / 2, 2 * (N2 - k2), n1 / 2, 2 * k2);
			}
		}
		if (N1 % 2 != 0) {
			const int n1 = N1 - 1;
			const int nkr = N2 * n1 + k2;
			const int nki = N2 * n1 + N2 - k2;
			if (k2 == 0 || (k2 == N2 / 2 && N2 % 2 == 0)) {
				print("x[%i] = z1[%i];\n", nkr, k2);
			} else if (k2 < N2 / 2 + 1) {
				print("x[%i] = z1[%i];\n", nkr, k2);
				print("x[%i] = z1[%i];\n", nki, N2 - k2);
			}
		}
	}
	for (int k2 = 0; k2 < N2 / 2 + 1; k2++) {
		for (int n1 = 0; n1 < N1; n1++) {
			const int kr = N2 * n1 + k2;
			const int ki = N2 * n1 + (N2 - k2);
			const auto W = twiddle(n1 * k2, N);
			if (k2 == 0) {

			} else if (k2 == N2 / 2 && N2 % 2 == 0) {

			} else {
				print("tmp0 = x[%i];\n", kr);
				print("x[%i] = x[%i] * %.17e + x[%i] * %.17e;\n", kr, kr, W.real(), ki, -W.imag());
				print("x[%i] = x[%i] * %.17e + tmp0 * %.17e;\n", ki, ki, W.real(), W.imag());
			}
		}
	}
	print("{\n");
	indent();
	for (int k1 = 0; k1 < N1; k1++) {
		const int kr = N2 * k1;
		print("const auto zr%i = x[%i];\n", k1, kr);
		print("constexpr double zi%i = 0;\n", k1);
	}
	switch (N1) {
	case 1:
		return;
	case 2:
		print("x[%i] = zr0 + zr1;\n", 0);
		print("x[%i] = zr0 - zr1;\n", N2);
		break;
	case 3:
		print("const auto tr1 = zr1 + zr2;\n");
		print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
		print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
		print("x[%i] = zr0 + tr1;\n", 0);
		print("x[%i] = tr2;\n", N2);
		print("x[%i] = -tr3;\n", N - N2);
		break;
	case 4:
		print("const auto tr1 = zr0 + zr2;\n");
		print("const auto tr2 = zr1 + zr3;\n");
		print("const auto tr3 = zr0 - zr2;\n");
		print("const auto tr4 = zr1 - zr3;\n");
		print("x[%i] = tr1 + tr2;\n", 0);
		print("x[%i] = tr3 ;\n", N2);
		print("x[%i] = - tr4;\n", N - N2);
		print("x[%i] = tr1 - tr2;\n", 2 * N2);
		print("x[%i] = tr3;\n", N - (3 * N2));
		print("x[%i] = -tr4;\n", 3 * N2);
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
		print("x[%i] = zr0 + tr5;\n", 0);
		print("x[%i] = tr8;\n", N2);
		print("x[%i] = -tr10;\n", N - N2);
		print("x[%i] = tr9;\n", 2 * N2);
		print("x[%i] = -tr11;\n", N - 2 * N2);
		break;
	case 6:
		print("const auto tr1 = zr2 + zr4;\n");
		print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
		print("const auto tr3 = (%24.17e) * (zr2 - zr4);\n", sin(M_PI / 3.0));
		print("const auto tr4 = zr5 + zr1;\n");
		print("const auto tr5 = std::fma(tr4, -0.5, zr3);\n");
		print("const auto tr6 = (%24.17e) * (zr5 - zr1);\n", sin(M_PI / 3.0));
		print("const auto tr7 = zr0 + tr1;\n");
		print("const auto tr8 = tr2;\n");
		print("const auto ti8 = -tr3;\n");
		print("const auto tr9 = tr2;\n");
		print("const auto ti9 = tr3;\n");
		print("const auto tr10 = zr3 + tr4;\n");
		print("const auto tr11 = tr5;\n");
		print("const auto ti11 = -tr6;\n");
		print("const auto tr12 = tr5;\n");
		print("const auto ti12 = tr6;\n");
		print("x[%i] = tr7 + tr10;\n", 0);
		print("x[%i] = tr8 - tr11;\n", N2);
		print("x[%i] = ti8 - ti11;\n", N - (N2));
		print("x[%i] = tr9 + tr12;\n", 2 * N2);
		print("x[%i] = ti9 + ti12;\n", N - (2 * N2));
		print("x[%i] = tr7 - tr10;\n", N - (3 * N2));

		break;

	default:
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			print("const auto txp%i = zr%i + zr%i;\n", j, j, N1 - j);
		}
		for (int j = 1; j <= (N1 - 1) / 2; j++) {
			print("const auto txm%i = zr%i - zr%i;\n", j, j, N1 - j);
		}
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			print("auto ap%i = zr0;\n", i);
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
				if (j == 1) {
					print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
				} else {
					print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
				}
			}
		}
		print("x[%i] = ", 0);
		for (int i = 0; i < N1; i++) {
			print_notab(" + zr%i", i);
		}
		print_notab(";\n");
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			print("x[%i] = ap%i;\n", i * N2, i);
			if (N - i * N2 != i * N2) {
				print("x[%i] = -bm%i;\n", N - (i * N2), i);
			}
			print("x[%i] = ap%i;\n", N - ((N1 - i) * N2), i);
			if (N - (N1 - i) * N2 != (N1 - i) * N2) {
				print("x[%i] = -bm%i;\n", (N1 - i) * N2, i);
			}
		}
	}
	deindent();
	print("}\n");

	for (int k2 = 1; k2 < (N2 + 1) / 2; k2++) {
		print("{\n");
		indent();
		for (int k1 = 0; k1 < N1; k1++) {
			const int kr = N2 * k1 + k2;
			const int ki = N2 * k1 + (N2 - k2);
			if (k2 == 0) {
				print("const auto zr%i = x[%i];\n", k1, kr);
				print("constexpr double zi%i = 0;\n", k1);
			} else if (k2 == N2 / 2 && N2 % 2 == 0) {
				const auto W = twiddle(k1 * k1, N1);
				print("const auto zr%i = x[%i];\n", k1, kr);
				print("constexpr double zi%i = 0;\n", k1);
			} else {
				const auto W = twiddle(k1 * k1, N1);
				print("const auto zr%i = x[%i];\n", k1, kr);
				print("const auto zi%i = x[%i];\n", k1, ki);
			}
		}

		switch (N1) {
		case 1:
			return;
		case 2:
			print("x[%i] = zr0 + zr1;\n", k2);
			print("x[%i] = zi0 + zi1;\n", N - k2);
			print("x[%i] = zr0 - zr1;\n", N - (N2 + k2));
			print("x[%i] = zi1 - zi0;\n", N2 + k2);
			break;
		case 3:
			print("const auto tr1 = zr1 + zr2;\n");
			print("const auto tr2 = std::fma(tr1, -0.5, zr0);\n");
			print("const auto tr3 = (%24.17e) * (zr1 - zr2);\n", sqrt(3) * 0.5);
			print("const auto ti1 = zi1 + zi2;\n");
			print("const auto ti2 = std::fma(ti1, -0.5, zi0);\n");
			print("const auto ti3 = (%24.17e) * (zi1 - zi2);\n", sqrt(3) * 0.5);
			print("x[%i] = zr0 + tr1;\n", k2);
			print("x[%i] = zi0 + ti1;\n", N - k2);
			print("x[%i] = tr2 + ti3;\n", N2 + k2);
			print("x[%i] = ti2 - tr3;\n", N - (N2 + k2));
			print("x[%i] = tr2 - ti3;\n", N - (2 * N2 + k2));
			print("x[%i] = -ti2 - tr3;\n", 2 * N2 + k2);
			break;
		case 4:
			print("const auto tr1 = zr0 + zr2;\n");
			print("const auto ti1 = zi0 + zi2;\n");
			print("const auto tr2 = zr1 + zr3;\n");
			print("const auto ti2 = zi1 + zi3;\n");
			print("const auto tr3 = zr0 - zr2;\n");
			print("const auto ti3 = zi0 - zi2;\n");
			print("const auto tr4 = zr1 - zr3;\n");
			print("const auto ti4 = zi1 - zi3;\n");
			print("x[%i] = tr1 + tr2;\n", k2);
			print("x[%i] = ti1 + ti2;\n", N - k2);
			print("x[%i] = tr3 + ti4;\n", N2 + k2);
			print("x[%i] = ti3 - tr4;\n", N - (N2 + k2));
			print("x[%i] = tr1 - tr2;\n", N - (2 * N2 + k2));
			print("x[%i] = ti2 - ti1 ;\n", (2 * N2 + k2));
			print("x[%i] = tr3 - ti4;\n", N - (3 * N2 + k2));
			print("x[%i] = -(ti3 + tr4);\n", 3 * N2 + k2);
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
			print("x[%i] = zr0 + tr5;\n", k2);
			print("x[%i] = zi0 + ti5;\n", N - k2);
			print("x[%i] = tr8 + ti10;\n", N2 + k2);
			print("x[%i] = ti8 - tr10;\n", N - (N2 + k2));
			print("x[%i] = tr9 + ti11;\n", 2 * N2 + k2);
			print("x[%i] = ti9 - tr11;\n", N - (2 * N2 + k2));
			print("x[%i] = tr9 - ti11;\n", N - (3 * N2 + k2));
			print("x[%i] = -ti9 - tr11;\n", (3 * N2 + k2));
			print("x[%i] = tr8 - ti10;\n", N - (4 * N2 + k2));
			print("x[%i] = -ti8 - tr10;\n", (4 * N2 + k2));
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
			print("x[%i] = tr7 + tr10;\n", 0);
			print("x[%i] = tr8 - tr11;\n", N2 + k2);
			print("x[%i] = tr9 + tr12;\n", 2 * N2 + k2);
			print("x[%i] = tr7 - tr10;\n", N - (3 * N2 + k2));
			print("x[%i] = tr8 + tr11;\n", N - (4 * N2 + k2));
			print("x[%i] = tr9 - tr12;\n", N - (5 * N2 + k2));
			print("x[%i] = ti7 + ti10;\n", N - k2);
			print("x[%i] = ti8 - ti11;\n", N - (N2 + k2));
			print("x[%i] = ti9 + ti12;\n", N - (2 * N2 + k2));
			print("x[%i] = -ti7 + ti10;\n", 3 * N2 + k2);
			print("x[%i] = -ti8 - ti11;\n", 4 * N2 + k2);
			print("x[%i] = -ti9 + ti12;\n", 5 * N2 + k2);
			break;
		default:
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("const auto txp%i = zr%i + zr%i;\n", j, j, N1 - j);
			}
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("const auto txm%i = zr%i - zr%i;\n", j, j, N1 - j);
			}
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("const auto typ%i = zi%i + zi%i;\n", j, j, N1 - j);
			}
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("const auto tym%i = zi%i - zi%i;\n", j, j, N1 - j);
			}
			for (int i = 1; i <= (N1 - 1) / 2; i++) {
				print("auto ap%i = zr0;\n", i);
				print("auto bp%i = zi0;\n", i);
				for (int j = 1; j <= (N1 - 1) / 2; j++) {
					print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
					print("bp%i = std::fma(typ%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
					if (j == 1) {
						print("double am%i = tym%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
						print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
					} else {
						print("am%i = std::fma(tym%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
						print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
					}
				}
			}
			print("x[%i] = ", k2);
			for (int i = 0; i < N1; i++) {
				print_notab(" + zr%i", i);
			}
			print_notab(";\n");
			print("x[%i] = ", N - k2);
			for (int i = 0; i < N1; i++) {
				print_notab(" + zi%i", i);
			}
			print_notab(";\n");
			for (int i = 1; i <= (N1 - 1) / 2; i++) {
				print("x[%i] = ap%i + am%i;\n", k2 + i * N2, i, i);
				print("x[%i] = bp%i - bm%i;\n", N - (k2 + i * N2), i, i);
				print("x[%i] = ap%i - am%i;\n", N - (k2 + (N1 - i) * N2), i, i);
				print("x[%i] = -bp%i - bm%i;\n", k2 + (N1 - i) * N2, i, i);
			}
		}
		deindent();
		print("}\n");

	}
	if (N2 % 2 == 0) {

		int k2 = N2 / 2;
		print("{\n");
		indent();
		for (int k1 = 0; k1 < N1; k1++) {
			const int kr = N2 * k1 + k2;
			print("const auto x%i = x[%i];\n", k1, kr);
		}

		switch (N1) {
		case 1:
			return;
		case 2:
			print("x[%i] = x0;\n", k2);
			print("x[%i] = -x1;\n", N - k2);
			break;
		case 3:
			print("const auto t1 = x1 - x2;\n");
			print("x[%i] = x0 + 0.5 * t1;\n", k2);
			print("x[%i] = x0 - t1;\n", N - k2);
			print("x[%i] = (%.17e) * (x1 + x2);\n", N2 + k2, -sin(M_PI / 3.0));
			break;
		case 4:
			print("const auto t1 = (x1 - x3) * (%.17e);\n", 1.0 / sqrt(2));
			print("const auto t2 = (x1 + x3) * (%.17e);\n", 1.0 / sqrt(2));
			print("x[%i] = x0 + t1;\n", k2);
			print("x[%i] = -x2 - t2;\n", N - k2);
			print("x[%i] = x0 - t1;\n", N2 + k2);
			print("x[%i] = x2 - t2;\n", N - (N2 + k2));
			break;
		case 5:
			print("const auto t1 = x1 - x4;\n");
			print("const auto t2 = x1 + x4;\n");
			print("const auto t3 = x2 - x3;\n");
			print("const auto t4 = x2 + x3;\n");
			print("const auto t5 = t1 - t3;\n");
			print("const auto t6 = x0 + 0.25 * t5;\n");
			print("const auto t7 = (%.17e) * (t1 + t3);\n", sqrt(5) * 0.25);
			print("x[%i] = t6 + t7;\n", k2);
			print("x[%i] = (%.17e) * t2 + (%.17e) * t4;\n", N - k2, -sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
			print("x[%i] = t6 - t7;\n", N2 + k2);
			print("x[%i] = (%.17e) * t2 + (%.17e) * t4;\n", N - (N2 + k2), -sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
			print("x[%i] = x0 - t5;\n", 2 * N2 + k2);
			break;
		case 6:
			print("const auto t1 = (%.17e) * (x5 - x1);\n", sin(M_PI / 3.0));
			print("const auto t2 = (%.17e) * (x2 + x4);\n", sin(M_PI / 3.0));
			print("const auto t3 = x2 - x4;\n");
			print("const auto t4 = x1 + x5;\n");
			print("const auto t5 = x0 + 0.5 * t3;\n");
			print("const auto t6 = -x3 - 0.5 * t4;\n");
			print("x[%i] = t5 - t1;\n", 0);
			print("x[%i] = t6 - t2;\n", N - k2);
			print("x[%i] = x0 - t3;\n", N2 + k2);
			print("x[%i] = x3 - t4;\n", N - (N2 + k2));
			print("x[%i] = t5 + t1;\n", 2 * N2 + k2);
			print("x[%i] = t6 + t2;\n", N - (2 * N2 + k2));
			break;
		}
		deindent();
		print("}\n");

	}

	deindent();
	print("}\n");

}
