#include "xorshift.hpp"
static std::random_device rd;

static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}


static const xorshift::state_type s_default_seed = {
  123456789, 362436069
};

xorshift::xorshift(void)
  : state_({rd(), rd()})
{}



void xorshift::seed(const state_type &seed) {
  state(seed);
}

void xorshift::seed(void) {
  state(s_default_seed);
}

const xorshift::state_type &xorshift::state(void) const {
  return state_;
}

void xorshift::state(const state_type &state) {
  state_ = state;
}

xorshift::result_type xorshift::min(void) {
  
  return std::numeric_limits<xorshift::result_type>::min();
}

xorshift::result_type xorshift::max(void) {
  return std::numeric_limits<xorshift::result_type>::max();
}


xorshift::result_type xorshift::operator()(void) {
  const xorshift::result_type s0 = state_.x;
  xorshift::result_type s1 = state_.y;
  const xorshift::result_type  result = s0 + s1;
  s1 ^= s0;
  state_.x = rotl(s0, 27) ^ s1 ^ (s1 << 7); // a, b
  state_.y = rotl(s1, 20); // c

  return result;
}

