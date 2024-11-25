// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

bool ieeeck(final int ISPEC, final double ZERO, final double ONE) {
  double NEGINF, NEGZRO, NEWZRO, POSINF;

  POSINF = ONE / ZERO;
  if (POSINF <= ONE) {
    return false;
  }

  NEGINF = -ONE / ZERO;
  if (NEGINF >= ZERO) {
    return false;
  }

  NEGZRO = ONE / (NEGINF + ONE);
  if (NEGZRO != ZERO) {
    return false;
  }

  NEGINF = ONE / NEGZRO;
  if (NEGINF >= ZERO) {
    return false;
  }

  NEWZRO = NEGZRO + ZERO;
  if (NEWZRO != ZERO) {
    return false;
  }

  POSINF = ONE / NEWZRO;
  if (POSINF <= ONE) {
    return false;
  }

  NEGINF *= POSINF;
  if (NEGINF >= ZERO) {
    return false;
  }

  POSINF *= POSINF;
  if (POSINF <= ONE) {
    return false;
  }

  // Return if we were only asked to check infinity arithmetic

  if (ISPEC == 0) return true;

  final NAN1 = POSINF + NEGINF;
  final NAN2 = POSINF / NEGINF;
  final NAN3 = POSINF / POSINF;
  final NAN4 = POSINF * ZERO;
  final NAN5 = NEGINF * NEGZRO;
  final NAN6 = NAN5 * ZERO;

  if (NAN1 == NAN1) {
    return false;
  }

  if (NAN2 == NAN2) {
    return false;
  }

  if (NAN3 == NAN3) {
    return false;
  }

  if (NAN4 == NAN4) {
    return false;
  }

  if (NAN5 == NAN5) {
    return false;
  }

  if (NAN6 == NAN6) {
    return false;
  }

  return true;
}
