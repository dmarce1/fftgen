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

void include(const char* name) {
	FILE* fp1 = fopen(name, "rt");
	if (fp1 == NULL) {
		printf("Fatal error could not include %s.\n", name);
		abort();
	}
	char buffer[1024];
	while (!feof(fp1)) {
		if( fgets(buffer, 1023, fp1) == 0 ) {
			break;
		}
		fprintf(fp, "%s", buffer);
	}

	fclose(fp1);

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
