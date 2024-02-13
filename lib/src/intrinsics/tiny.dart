import 'dart:typed_data';

double tiny(final double _) {
  // minimum **normal** positive 64 bits float (exp=1, m=0)
  return Uint32List.fromList([0x00000000, 0x00100000])
      .buffer
      .asFloat64List(0)
      .first;
}
