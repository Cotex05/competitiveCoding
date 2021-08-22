int hamming(string a, string b) {
  int d = 0;
  for (int i = 0; i < k; i++) {
    if (a[i] != b[i]) d++;
  }
    return d;
}


int hamming(int a, int b) {
  return __builtin_popcount(a^b);
}
