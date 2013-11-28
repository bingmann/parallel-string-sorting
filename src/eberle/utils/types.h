#ifndef TYPES_H_
#define TYPES_H_

namespace types {

typedef unsigned char* string;
typedef unsigned int UINT;

struct AS {
	UINT lcp;
	string text;

	AS() {
		lcp = 0;
		text = 0;
	}

	AS(UINT lcp, string text) {
		this->lcp = lcp;
		this->text = text;
	}
};

} // namespace types

#endif // TYPES_H_
