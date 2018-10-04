#include "shared.h"

void Assert(bool condition, const char * s)
{
	if (!condition) {
		printf("assert failed: %s", s);
		throw s;
	}
}
