#ifndef HOOKS_H
#define HOOKS_H

#include "prism/domain.h"
#include "speclib/speclib.h"

extern BSystem* (*user_build) (Field*,Bedge*,const char*,double);

extern void (*user_advect)(Domain*);

extern void (*user_transform)();
extern void (*user_gradz)();

#endif

