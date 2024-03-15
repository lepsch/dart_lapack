extension RangeExtension on int {
  /// Returns a generator from [this] up to, but not including, the [end].
  Iterable<int> to(int end, {int step = 1}) sync* {
    for (var i = this; step < 0 ? i > end : i < end; i += step) {
      yield i;
    }
  }

  /// Returns a generator from [this] through the [end].
  Iterable<int> through(int end, {int step = 1}) sync* {
    for (var i = this; step < 0 ? i >= end : i <= end; i += step) {
      yield i;
    }
  }
}
