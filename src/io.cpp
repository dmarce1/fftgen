#include "fftgen.hpp"


FILE* fp = nullptr;

int tabs = 0;

void indent() {
	tabs++;
}

void deindent() {
	tabs--;
}

FILE* get_fp() {
	return fp;
}

int get_ntabs() {
	return tabs;
}

void set_file(std::string str) {
	if (fp) {
		fclose(fp);
	}
	if (str.empty()) {
		fp = nullptr;
	} else {
		str = "./gen/" + str;
		fp = fopen(str.c_str(), "wt");
	}
}
