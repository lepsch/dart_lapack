      bool    FUNCTION LSAME(CA,CB);

*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String    CA,CB;
      // ..

* =====================================================================

      // .. Intrinsic Functions ..
      // INTRINSIC ICHAR
      // ..
      // .. Local Scalars ..
      int     INTA,INTB,ZCODE;
      // ..

      // Test if the characters are equal

      LSAME = CA == CB
      if (LSAME) RETURN;

      // Now test for equivalence if both characters are alphabetic.

      ZCODE = ICHAR('Z')

      // Use 'Z' rather than 'A' so that ASCII can be detected on Prime
      // machines, on which ICHAR returns a value with bit 8 set.
      // ICHAR('A') on Prime machines returns 193 which is the same as
      // ICHAR('A') on an EBCDIC machine.

      INTA = ICHAR(CA)
      INTB = ICHAR(CB)

      if (ZCODE == 90 || ZCODE == 122) {

         // ASCII is assumed - ZCODE is the ASCII code of either lower or
         // upper case 'Z'.

          if (INTA.GE.97 && INTA.LE.122) INTA = INTA - 32;
          if (INTB.GE.97 && INTB.LE.122) INTB = INTB - 32;

      } else if (ZCODE == 233 || ZCODE == 169) {

         // EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
         // upper case 'Z'.

          if (INTA.GE.129 && INTA.LE.137 || INTA.GE.145 && INTA.LE.153 || INTA.GE.162 && INTA.LE.169) INTA = INTA + 64           IF (INTB.GE.129 && INTB.LE.137 || INTB.GE.145 && INTB.LE.153 || INTB.GE.162 && INTB.LE.169) INTB = INTB + 64;

      } else if (ZCODE == 218 || ZCODE == 250) {

         // ASCII is assumed, on Prime machines - ZCODE is the ASCII code
         // plus 128 of either lower or upper case 'Z'.

          if (INTA.GE.225 && INTA.LE.250) INTA = INTA - 32;
          if (INTB.GE.225 && INTB.LE.250) INTB = INTB - 32;
      }
      LSAME = INTA == INTB

      // RETURN

      // End of LSAME

      }
