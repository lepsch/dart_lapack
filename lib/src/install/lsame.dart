      bool lsame(CA, CB ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             CA, CB;
      // ..

// =====================================================================

      // .. Intrinsic Functions ..
      // INTRINSIC ICHAR
      int                INTA, INTB, ZCODE;

      // Test if the characters are equal

      lsame = CA == CB;
      if (lsame) return;

      // Now test for equivalence if both characters are alphabetic.

      ZCODE = ICHAR( 'Z' );

      // Use 'Z' rather than 'A' so that ASCII can be detected on Prime
      // machines, on which ICHAR returns a value with bit 8 set.
      // ICHAR('A') on Prime machines returns 193 which is the same as
      // ICHAR('A') on an EBCDIC machine.

      INTA = ICHAR( CA );
      INTB = ICHAR( CB );

      if ( ZCODE == 90 || ZCODE == 122 ) {

         // ASCII is assumed - ZCODE is the ASCII code of either lower or
         // upper case 'Z'.

         if (INTA >= 97 && INTA <= 122) INTA = INTA - 32;
         if (INTB >= 97 && INTB <= 122) INTB = INTB - 32;

      } else if ( ZCODE == 233 || ZCODE == 169 ) {

         // EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
         // upper case 'Z'.

         if (INTA >= 129 && INTA <= 137 || INTA >= 145 && INTA <= 153 || INTA >= 162 && INTA <= 169) INTA = INTA + 64;
         IF( INTB >= 129 && INTB <= 137 || INTB >= 145 && INTB <= 153 || INTB >= 162 && INTB <= 169 ) INTB = INTB + 64;

      } else if ( ZCODE == 218 || ZCODE == 250 ) {

         // ASCII is assumed, on Prime machines - ZCODE is the ASCII code
         // plus 128 of either lower or upper case 'Z'.

         if (INTA >= 225 && INTA <= 250) INTA = INTA - 32;
         if (INTB >= 225 && INTB <= 250) INTB = INTB - 32;
      }
      lsame = INTA == INTB;

      // }
