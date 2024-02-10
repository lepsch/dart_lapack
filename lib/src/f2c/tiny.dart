import 'dart:typed_data';

double tiny(final double _) {
  // minimum **normal** positive double (exp=1, m=0)
  return Uint64List.fromList([0x10000000000000])
      .buffer
      .asFloat64List(0)
      .first;
}
