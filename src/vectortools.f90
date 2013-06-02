module vectortools
      implicit none
      contains
        function vectorabs(x)
		real (kind=4):: x(3), vectorabs
		vectorabs = sqrt(x(1)**2. + x(2)**2. +x(3)**2.)
        end function vectorabs
        function crossproduct(x, y)
		real (kind=4)::  crossproduct(1:3), x(1:3), y(1:3)
		crossproduct(1) = x(2)*y(3) - x(3)*y(2)
		crossproduct(2) = x(3)*y(1) - x(1)*y(3)
		crossproduct(3) = x(1)*y(2) - x(2)*y(1)
	end function crossproduct
		
	function dotproduct(x, y)
		real (kind=4) ::  dotproduct, x(1:3), y(1:3)
		integer :: i
		dotproduct = x(1)*y(1)+x(2)*y(2)+x(3)*y(3)	
	end function dotproduct

	function finddistance(x_1,x_2,point)
		real (kind=4):: finddistance
		real (kind=4):: point(1:3), x_1(1:3), x_2(1:3)
		real (kind=4):: dummy1(1:3), dummy2(1:3)
		dummy1 = crossproduct((point - x_1), (point -x_2) )
		dummy2 = x_2 - x_1
		finddistance = vectorabs(dummy1) / vectorabs(dummy2)
	end function finddistance
	
	function linedistance(init_point,target_point,halo_centre)
		real (kind=4) :: init_point(3),target_point(3),halo_centre(3)
		real (kind=4) :: linedistance
		linedistance = dotproduct((target_point-init_point),(halo_centre-init_point)) / vectorabs(target_point-init_point)		
	end function linedistance
end module vectortools
