#include "fftgen.hpp"

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
	case 15:

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
/*
 void print_special_short_fft(int r, std::vector<std::string> in_, std::vector<std::string> out_) {
 std::vector<const char*> in(in_.size());
 std::vector<const char*> out(out_.size());
 for (int i = 0; i < in.size(); i++) {
 in[i] = in_[i].c_str();
 out[i] = out_[i].c_str();
 }
 switch (r) {
 case 2:
 print("%s = %s;\n", out[0], in[0]);
 break;
 case 3:
 print("const auto& tr1 = %s;\n", in[1]);
 print("const auto tr2 = std::fma(tr1, -0.5, %s);\n", in[0]);
 print("const auto tr3 = (%24.17e) * %s;\n", sqrt(3) * 0.5, in[1]);
 print("%s = tr2;\n", out[0]);
 print("%s = -tr3;\n", out[1]);
 break;
 case 4:
 print("const auto& tr1 = %s;\n", in[0]);
 print("const auto& tr3 = %s;\n", in[0]);
 print("const auto& tr2 = %s;\n", in[1]);
 print("%s = tr3;\n", out[0]);
 print("%s = tr1 - tr2;\n", out[1]);
 break;
 case 5:
 print("const auto tr1 = %s;\n", in[1]);
 print("const auto tr2 = %s;\n", in[2]);
 print("const auto tr3 = %s;\n", in[1]);
 print("const auto tr4 = %s;\n", in[2]);
 print("const auto tr6 = (%24.17e) * (tr1 - tr2);\n", sqrt(5) * 0.25);
 print("const auto tr7 = std::fma(tr5, -0.25, %s);\n", in[0]);
 print("const auto tr8 = tr7 + tr6;\n");
 print("const auto tr9 = tr7 - tr6;\n");
 print("const auto tr10 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 5.0), sin(2.0 * M_PI / 10.0));
 print("const auto tr11 = std::fma((%24.17e), tr3, (%24.17e) * tr4);\n", sin(2.0 * M_PI / 10.0), -sin(2.0 * M_PI / 5.0));
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
 case 15:

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
 }*/

void print_skew_short_fft(int N1, std::vector<std::string> in_, std::vector<std::string> out_) {
	std::vector<const char*> in(in_.size());
	std::vector<const char*> out(out_.size());
	for (int i = 0; i < in.size(); i++) {
		in[i] = in_[i].c_str();
		out[i] = out_[i].c_str();
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
	case 7: {
		for (int k1 = 0; k1 < N1; k1++) {
			print("const auto x%i = %s;\n", k1, in[k1]);
		}
		constexpr double theta = 2.0 * M_PI / 7.0;
		const double c1 = cos(theta);
		const double c2 = cos(2.0 * theta);
		const double c3 = cos(3.0 * theta);
		const double c4 = sin(theta);
		const double c5 = sin(2.0 * theta);
		const double c6 = sin(3.0 * theta);
		print("const auto t1 = x0 - 0.5 * x6;\n");
		print("const auto t2 = x4 - x6;\n");
		print("const auto t3 = x2 - x6;\n");
		print("const auto t4 = -0.5 * x1;\n");
		print("const auto t5 = x5 - x1;\n");
		print("const auto t6 = x3 - x1;\n");
		print("const auto er0 = x0 + x2 + x4 + x6;\n");
		print("const auto er1 = std::fma(t3, %.17e, std::fma(t2, %.17e, t1));\n", c1, c2);
		print("const auto ei1 = std::fma(x2, %.17e, std::fma(x4, %.17e, x6 * %.17e));\n", -c4, -c5, -c6);
		print("const auto er3 = std::fma(t3, %.17e, std::fma(t2, %.17e, t1));\n", c3, c1);
		print("const auto ei3 = std::fma(x2, %.17e, std::fma(x4, %.17e, x6 * %.17e));\n", -c6, c4, -c5);
		print("const auto er5 = std::fma(t3, %.17e, std::fma(t2, %.17e, t1));\n", c2, c3);
		print("const auto ei5 = std::fma(x2, %.17e, std::fma(x4, %.17e, x6 * %.17e));\n", c5, -c6, -c4);
		print("const auto or0 = x5 + x3 + x1;\n");
		print("const auto or1 = std::fma(t5, %.17e, std::fma(t6, %.17e, t4));\n", c1, c2);
		print("const auto oi1 = std::fma(x5, %.17e, std::fma(x3, %.17e, x1 * %.17e));\n", c4, c5, c6);
		print("const auto or3 = std::fma(t5, %.17e, std::fma(t6, %.17e, t4));\n", c3, c1);
		print("const auto oi3 = std::fma(x5, %.17e, std::fma(x3, %.17e, x1 * %.17e));\n", c6, -c4, c5);
		print("const auto or5 = std::fma(t5, %.17e, std::fma(t6, %.17e, t4));\n", c2, c3);
		print("const auto oi5 = std::fma(x5, %.17e, std::fma(x3, %.17e, x1 * %.17e));\n", -c5, c6, c4);
		print("%s = er1 - or1;\n", out[0]);
		print("%s = ei1 - oi1;\n", out[6]);
		print("%s = er3 - or3;\n", out[1]);
		print("%s = ei3 - oi3;\n", out[5]);
		print("%s = er5 - or5;\n", out[2]);
		print("%s = ei5 - oi5;\n", out[4]);
		print("%s = er0 - or0;\n", out[3]);

	}
		break;
	case 8: {
		for (int k1 = 0; k1 < N1; k1++) {
			print("const auto x%i = %s;\n", k1, in[k1]);
		}
		const auto tw1 = twiddle(1, 16);
		const auto tw3 = twiddle(3, 16);
		const auto tw5 = twiddle(5, 16);
		const auto tw7 = twiddle(7, 16);
		print("const auto tr1 = (%.17e) * (x2 - x6);\n", 1.0 / sqrt(2));
		print("const auto tr2 = (%.17e) * (x2 + x6);\n", 1.0 / sqrt(2));
		print("const auto tr3 = (%.17e) * (x3 - x7);\n", 1.0 / sqrt(2));
		print("const auto tr4 = (%.17e) * (x3 + x7);\n", 1.0 / sqrt(2));
		print("const auto or1 = x1 + tr3;\n");
		print("const auto or3 = x1 - tr3;\n");
		print("const auto oi3 = x5 - tr4;\n");
		print("const auto oi1 = -tr4 - x5;\n");
		print("const auto er1 = x0 + tr1;\n");
		print("const auto er3 = x0 - tr1;\n");
		print("const auto ei3 = x4 - tr2;\n");
		print("const auto ei1 = -tr2 - x4;\n");
		print("const auto or19 = (9.23879532511286738e-01) * or1;\n");
		print("const auto or39 = (9.23879532511286738e-01) * or3;\n");
		print("const auto oi19 = (9.23879532511286738e-01) * oi1;\n");
		print("const auto oi39 = (9.23879532511286738e-01) * oi3;\n");
		print("const auto or13 = (3.82683432365089837e-01) * or1;\n");
		print("const auto or33 = (3.82683432365089837e-01) * or3;\n");
		print("const auto oi13 = (3.82683432365089837e-01) * oi1;\n");
		print("const auto oi33 = (3.82683432365089837e-01) * oi3;\n");
		print("%s =  er1 + or19 + oi13;\n", out[0]);
		print("%s =  ei1 - or13 + oi19;\n", out[7]);
		print("%s =  er3 + or33 + oi39;\n", out[1]);
		print("%s =  ei3 - or39 + oi33;\n", out[6]);
		print("%s =  er3 - or33 - oi39;\n", out[2]);
		print("%s = -ei3 - or39 + oi33;\n", out[5]);
		print("%s =  er1 - or19 - oi13;\n", out[3]);
		print("%s = -ei1 - or13 + oi19;\n", out[4]);
	}
		break;
	case 9: {
		for (int k1 = 0; k1 < N1; k1++) {
			print("const auto x%i = %s;\n", k1, in[k1]);
		}
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
		print("const auto sr9 = x2 + x4;\n");
		print("const auto sr10 = sr9 + x8;\n");
		print("const auto sr11 = sr10 + x6;\n");
		print("const auto sr12 = sr11 + x0;\n");
		print("const auto sr13 = x2 - x4;\n");
		print("const auto sr14 = sr13 + x8;\n");
		print("const auto sr16 = x4 - x8;\n");
		print("const auto sr17 = x8 - x2;\n");
		print("const auto sr19 = -x4 - x8;\n");
		print("const auto mr1 = sr10 * (%.17e);\n", ar1);
		print("const auto mr3 = x6 * (%.17e);\n", ar3);
		print("const auto mi4 = x6 * (%.17e);\n", ai4);
		print("const auto mr5 = sr13 * (%.17e);\n", ar5);
		print("const auto mr6 = sr16 * (%.17e);\n", ar6);
		print("const auto mr7 = sr17 * (%.17e);\n", ar7);
		print("const auto mi8 = sr9 * (%.17e);\n", ai8);
		print("const auto mi9 = sr19 * (%.17e);\n", ai9);
		print("const auto mi10 = sr17 * (%.17e);\n", ai10);
		print("const auto sr21 = mr1 + mr1;\n");
		print("const auto sr22 = sr21 + mr1;\n");
		print("const auto sr26 = sr12 + mr3;\n");
		print("const auto sr27 = sr26 + sr21;\n");
		print("const auto sr28 = sr27 + mr5;\n");
		print("const auto sr30 = sr27 - mr6;\n");
		print("const auto sr32 = sr27 - mr5;\n");
		print("const auto si34 = mi4 + mi8;\n");
		print("const auto si36 = mi4 - mi9;\n");
		print("const auto si38 = mi4 - mi8;\n");
		print("const auto tr9 = x7 + x5;\n");
		print("const auto tr10 = tr9 + x1;\n");
		print("const auto tr11 = tr10 + x3;\n");
		print("const auto tr15 = x7 - x5;\n");
		print("const auto tr14 = -tr15 - x1;\n");
		print("const auto tr16 = x5 - x1;\n");
		print("const auto tr17 = x1 - x7;\n");
		print("const auto tr19 = x5 + x1;\n");
		print("const auto tr20 = -x1 + x7;\n");
		print("const auto nr1 = tr10 * (%.17e);\n", ar1);
		print("const auto nr3 = x3 * (%.17e);\n", ar3);
		print("const auto ni4 = x3 * (%.17e);\n", -ai4);
		print("const auto nr5 = tr15 * (%.17e);\n", ar5);
		print("const auto nr6 = tr16 * (%.17e);\n", ar6);
		print("const auto nr7 = tr17 * (%.17e);\n", ar7);
		print("const auto ni8 = tr9 * (%.17e);\n", -ai8);
		print("const auto ni9 = tr19 * (%.17e);\n", ai9);
		print("const auto ni10 = tr20 * (%.17e);\n", ai10);
		print("const auto tr21 = nr1 + nr1;\n");
		print("const auto tr22 = tr21 + nr1;\n");
		print("const auto tr26 = tr11 + nr3;\n");
		print("const auto tr27 = tr26 + tr21;\n");
		print("const auto tr28 = tr27 + nr5;\n");
		print("const auto tr30 = tr27 - nr6;\n");
		print("const auto tr32 = tr27 - nr5;\n");
		print("const auto ti34 = ni4 + ni8;\n");
		print("const auto ti36 = ni4 - ni9;\n");
		print("const auto ti38 = ni4 - ni8;\n");
		print("const auto er0 = sr12;\n");
		print("const auto ei3 = sr14 * (%.17e);\n", ai2);
		print("const auto er3 = sr22 + sr12;\n");
		print("const auto er1 = sr28 + mr6;\n");
		print("const auto er7 = sr30 + mr7;\n");
		print("const auto er5 = sr32 - mr7;\n");
		print("const auto ei1 = si34 + mi9;\n");
		print("const auto ei7 = si36 + mi10;\n");
		print("const auto ei5 = mi10 - si38;\n");
		print("const auto or0 = tr11;\n");
		print("const auto oi3 = tr14 * (%.17e);\n", ai2);
		print("const auto or3 = tr22 + tr11;\n");
		print("const auto or1 = tr28 + nr6;\n");
		print("const auto or7 = tr30 + nr7;\n");
		print("const auto or5 = tr32 - nr7;\n");
		print("const auto oi1 = ti34 + ni9;\n");
		print("const auto oi7 = ti36 + ni10;\n");
		print("const auto oi5 = ni10 - ti38;\n");
		print("%s = er1 - or1;\n", out[0]);
		print("%s = ei1 - oi1;\n", out[8]);
		print("%s = er3 - or3;\n", out[1]);
		print("%s = ei3 - oi3;\n", out[7]);
		print("%s = er5 - or5;\n", out[2]);
		print("%s = ei5 - oi5;\n", out[6]);
		print("%s = er7 - or7;\n", out[3]);
		print("%s = ei7 - oi7;\n", out[5]);
		print("%s = er0 - or0;\n", out[4]);

	}
		break;
	default:
		print("const auto zr0 = x[%i];\n", in[0]);
		for (int n1 = 1; n1 < N1; n1++) {
			const auto W = twiddle(n1, 2 * N1);
			print("const auto zr%i = %s * (%.17e);\n", n1, in[n1], W.real());
			print("const auto zi%i = %s * (%.17e);\n", n1, in[n1], W.imag());
		}
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
			for (int j = 1; j <= (N1 - 1) / 2; j++) {
				print("ap%i = std::fma(txp%i, (%24.17e), ap%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
				if (j == 1) {
					print("double bp%i = typ%i * (%24.17e);\n", i, j, cos(2.0 * M_PI * j * i / N1));
					print("double am%i = tym%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
					print("double bm%i = txm%i * (%24.17e);\n", i, j, sin(2.0 * M_PI * j * i / N1));
				} else {
					print("bp%i = std::fma(typ%i, (%24.17e), bp%i);\n", i, j, cos(2.0 * M_PI * j * i / N1), i);
					print("am%i = std::fma(tym%i, (%24.17e), am%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
					print("bm%i = std::fma(txm%i, (%24.17e), bm%i);\n", i, j, sin(2.0 * M_PI * j * i / N1), i);
				}
			}
		}
		print("%s = ", out[0]);
		for (int i = 0; i < N1; i++) {
			print_notab(" + zr%i", i);
		}
		print_notab(";\n");
		print("%s = ", out[N1 - 1]);
		for (int i = 1; i < N1; i++) {
			print_notab(" + zi%i", i);
		}
		print_notab(";\n");
		for (int i = 1; i <= (N1 - 1) / 2; i++) {
			if (i == (N1 - 1) / 2) {
				print("%s = ap%i + am%i;\n", out[2 * i], i, i);
			} else {
				print("%s = ap%i + am%i;\n", out[2 * i], i, i);
				print("%s = bp%i - bm%i;\n", out[2 * i + 1], i, i);
			}
		}
		break;
	}
}

