module mpas_typecast

   use mpas_kind_types
   use mpas_abort, only : mpas_dmpar_global_abort

   implicit none

   interface encode
      module procedure encode_character
      module procedure encode_ikind
      module procedure encode_i8kind
      module procedure encode_r4kind
      module procedure encode_r8kind
   end interface encode

   interface decode
      module procedure decode_character
      module procedure decode_ikind
      module procedure decode_i8kind
      module procedure decode_r4kind
      module procedure decode_r8kind
   end interface decode

contains

   subroutine check_buffersize(buffersize, target)

      implicit none

      integer, intent(in) :: buffersize, target

      if (target > buffersize) then
         call mpas_dmpar_global_abort("Buffer overflow while encoding/decoding character buffer")
      end if
      !DIR$ ATTRIBUTES INLINE :: check_buffersize
   end subroutine check_buffersize

   subroutine encode_character(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(inout) :: buffer
      integer, intent(inout) :: cursor
      character(len=*), intent(in) :: value

      integer :: length, i

      length = len(value)
      call check_buffersize(size(buffer), cursor+length)
      do i=1,length
         buffer(cursor+i-1) = value(i:i)
      end do
      cursor = cursor + length
      !DIR$ ATTRIBUTES INLINE :: encode_character
   end subroutine encode_character

   subroutine encode_ikind(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(inout) :: buffer
      integer, intent(inout) :: cursor
      integer, intent(in) :: value

      type :: data_t
      integer :: x
      end type data_t
      type(data_t), target :: d
      character(len=1), dimension(:), allocatable :: enc
      character(len=1), dimension(1) :: encdummy
      integer :: length

      d % x = value
      length = size(transfer(d, encdummy))
      call check_buffersize(size(buffer), cursor+length)
      allocate(enc(length))
      enc = transfer(d, enc)
      buffer(cursor:cursor+length-1) = enc
      cursor = cursor + length
      deallocate(enc)
      !DIR$ ATTRIBUTES INLINE :: encode_ikind
   end subroutine encode_ikind

   subroutine encode_i8kind(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(inout) :: buffer
      integer, intent(inout) :: cursor
      integer (kind=I8KIND), intent(in) :: value

      type :: data_t
         integer (kind=I8KIND) :: x
      end type data_t
      type(data_t), target :: d
      character(len=1), dimension(:), allocatable :: enc
      character(len=1), dimension(1) :: encdummy
      integer :: length

      d % x = value
      length = size(transfer(d, encdummy))
      call check_buffersize(size(buffer), cursor+length)
      allocate(enc(length))
      enc = transfer(d, enc)
      buffer(cursor:cursor+length-1) = enc
      cursor = cursor + length
      deallocate(enc)
      !DIR$ ATTRIBUTES INLINE :: encode_i8kind
   end subroutine encode_i8kind

   subroutine encode_r4kind(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(inout) :: buffer
      integer, intent(inout) :: cursor
      real(KIND=R4KIND), intent(in) :: value

      type :: data_t
         real(KIND=R4KIND) :: x
      end type data_t
      type(data_t), target :: d
      character(len=1), dimension(:), allocatable :: enc
      character(len=1), dimension(1) :: encdummy
      integer :: length

      d % x = value
      length = size(transfer(d, encdummy))
      call check_buffersize(size(buffer), cursor+length)
      allocate(enc(length))
      enc = transfer(d, enc)
      buffer(cursor:cursor+length-1) = enc
      cursor = cursor + length
      deallocate(enc)
      !DIR$ ATTRIBUTES INLINE :: encode_r4kind
   end subroutine encode_r4kind

   subroutine encode_r8kind(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(inout) :: buffer
      integer, intent(inout) :: cursor
      real(KIND=R8KIND), intent(in) :: value

      type :: data_t
         real(KIND=R8KIND) :: x
      end type data_t
      type(data_t), target :: d
      character(len=1), dimension(:), allocatable :: enc
      character(len=1), dimension(1) :: encdummy
      integer :: length

      d % x = value
      length = size(transfer(d, encdummy))
      call check_buffersize(size(buffer), cursor+length)
      allocate(enc(length))
      enc = transfer(d, enc)
      buffer(cursor:cursor+length-1) = enc
      cursor = cursor + length
      deallocate(enc)
      !DIR$ ATTRIBUTES INLINE :: encode_r8kind
   end subroutine encode_r8kind

   subroutine decode_character(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(in) :: buffer
      integer, intent(inout) :: cursor
      character(len=*), intent(out) :: value

      integer :: length, i

      length = len(value)
      call check_buffersize(size(buffer), cursor+length)
      do i=1,length
         value(i:i) = buffer(cursor+i-1)
      end do
      cursor = cursor + length
      !DIR$ ATTRIBUTES INLINE :: decode_character
   end subroutine decode_character

   subroutine decode_ikind(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(in) :: buffer
      integer, intent(inout) :: cursor
      integer, intent(out) :: value

      type :: data_t
         integer :: x
      end type data_t
      type(data_t), target :: d
      character(len=1), dimension(:), allocatable :: enc
      character(len=1), dimension(1) :: encdummy
      integer :: length

      length = size(transfer(d, encdummy))
      call check_buffersize(size(buffer), cursor+length)
      allocate(enc(length))
      enc = buffer(cursor:cursor+length-1)
      d = transfer(enc, d)
      value = d % x
      cursor = cursor + length
      deallocate(enc)
      !DIR$ ATTRIBUTES INLINE :: decode_ikind
   end subroutine decode_ikind

   subroutine decode_i8kind(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(in) :: buffer
      integer, intent(inout) :: cursor
      integer (kind=I8KIND), intent(out) :: value

      type :: data_t
         integer (kind=I8KIND) :: x
      end type data_t
      type(data_t), target :: d
      character(len=1), dimension(:), allocatable :: enc
      character(len=1), dimension(1) :: encdummy
      integer :: length

      length = size(transfer(d, encdummy))
      call check_buffersize(size(buffer), cursor+length)
      allocate(enc(length))
      enc = buffer(cursor:cursor+length-1)
      d = transfer(enc, d)
      value = d % x
      cursor = cursor + length
      deallocate(enc)
      !DIR$ ATTRIBUTES INLINE :: decode_i8kind
   end subroutine decode_i8kind

   subroutine decode_r4kind(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(in) :: buffer
      integer, intent(inout) :: cursor
      real(KIND=R4KIND), intent(out) :: value

      type :: data_t
         real(KIND=R4KIND) :: x
      end type data_t
      type(data_t), target :: d
      character(len=1), dimension(:), allocatable :: enc
      character(len=1), dimension(1) :: encdummy
      integer :: length

      length = size(transfer(d, encdummy))
      call check_buffersize(size(buffer), cursor+length)
      allocate(enc(length))
      enc = buffer(cursor:cursor+length-1)
      d = transfer(enc, d)
      value = d % x
      cursor = cursor + length
      deallocate(enc)
      !DIR$ ATTRIBUTES INLINE :: decode_r4kind
   end subroutine decode_r4kind

   subroutine decode_r8kind(buffer, cursor, value)

      implicit none

      character(len=1), dimension(:), intent(in) :: buffer
      integer, intent(inout) :: cursor
      real(KIND=R8KIND), intent(out) :: value

      type :: data_t
         real(KIND=R8KIND) :: x
      end type data_t
      type(data_t), target :: d
      character(len=1), dimension(:), allocatable :: enc
      character(len=1), dimension(1) :: encdummy
      integer :: length

      length = size(transfer(d, encdummy))
      call check_buffersize(size(buffer), cursor+length)
      allocate(enc(length))
      enc = buffer(cursor:cursor+length-1)
      d = transfer(enc, d)
      value = d % x
      cursor = cursor + length
      deallocate(enc)
      !DIR$ ATTRIBUTES INLINE :: decode_r8kind
   end subroutine decode_r8kind

end module mpas_typecast
