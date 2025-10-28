subroutine dofler_mark(elem_err, mvolele, theta, mark)
  implicit none
  integer, intent(in) :: mvolele
  real(kind=8), intent(in) :: elem_err(mvolele)
  real(kind=8), intent(in) :: theta
  integer, intent(out) :: mark(mvolele)

  ! local arrays for sorting
  integer :: i
  real(kind=8), allocatable :: err_copy(:)
  integer, allocatable :: idx(:)
  real(kind=8) :: total_err, accum
  integer :: j, sel

  allocate(err_copy(mvolele))
  allocate(idx(mvolele))

  ! init
  total_err = 0.0_8
  do i = 1, mvolele
     err_copy(i) = elem_err(i)
     total_err = total_err + err_copy(i)
     mark(i) = 0
     idx(i) = i
  end do

  if (total_err <= 0.0_8) then
     ! nothing to refine
     return
  end if

  ! Simple selection: sort elements by error descending
  ! We implement a simple selection sort for clarity (mvolele can be large; you may replace with quicksort)
  integer :: p, tmpi
  real(kind=8) :: tmpr
  do i = 1, mvolele-1
     p = i
     do j = i+1, mvolele
        if (err_copy(j) > err_copy(p)) p = j
     end do
     if (p /= i) then
        tmpr = err_copy(i); err_copy(i) = err_copy(p); err_copy(p) = tmpr
        tmpi = idx(i); idx(i) = idx(p); idx(p) = tmpi
     end if
  end do

  accum = 0.0_8
  sel = 0
  do i = 1, mvolele
     accum = accum + err_copy(i)
     mark(idx(i)) = 1
     sel = sel + 1
     if (accum >= theta * total_err) exit
  end do

  ! now mark array has 1 for chosen elements
  deallocate(err_copy)
  deallocate(idx)
  return
end subroutine dofler_mark