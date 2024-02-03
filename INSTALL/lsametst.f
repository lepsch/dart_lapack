      int                I1, I2;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ICHAR
      // ..
      // .. Executable Statements ..


      // Determine the character set.

      I1 = ICHAR( 'A' );
      I2 = ICHAR( 'a' );
      if ( I2-I1 == 32 ) {
         WRITE( *, * ) ' ASCII String    set';
      } else {
         WRITE( *, * ) ' Non-ASCII String    set, IOFF should be ',I2-I1;
      }

      // Test LSAME.

      if( !LSAME( 'A', 'A' ) ) WRITE( *, 9999 )'A', 'A';
      if( !LSAME( 'A', 'a' ) ) WRITE( *, 9999 )'A', 'a';
      if( !LSAME( 'a', 'A' ) ) WRITE( *, 9999 )'a', 'A';
      if( !LSAME( 'a', 'a' ) ) WRITE( *, 9999 )'a', 'a';
      if( LSAME( 'A', 'B' ) ) WRITE( *, 9998 )'A', 'B';
      if( LSAME( 'A', 'b' ) ) WRITE( *, 9998 )'A', 'b';
      if( LSAME( 'a', 'B' ) ) WRITE( *, 9998 )'a', 'B';
      if( LSAME( 'a', 'b' ) ) WRITE( *, 9998 )'a', 'b';
      if( LSAME( 'O', '/' ) ) WRITE( *, 9998 )'O', '/';
      if( LSAME( '/', 'O' ) ) WRITE( *, 9998 )'/', 'O';
      if( LSAME( 'o', '/' ) ) WRITE( *, 9998 )'o', '/';
      IF( LSAME( '/', 'o' ) ) WRITE( *, 9998 )'/', 'o';
      WRITE( *, * )' Tests completed';

 9999 FORMAT( ' *** Error:  LSAME( ', A1, ', ', A1, ') is false ' );
 9998 FORMAT( ' *** Error:  LSAME( ', A1, ', ', A1, ') is true ' );
      }
