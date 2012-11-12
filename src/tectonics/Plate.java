package tectonics;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Plate
{
	private static boolean DEBUG = true; //TODO: make final
	
	private static final int INITIAL_SPEED_X = 1;
	private static final int DEFORMATION_WEIGHT = 5;
	private static final float CONT_BASE = 1.0f; ///< Height limit that separates seas from dry land.
	
	/*
	// http://en.wikipedia.org/wiki/Methods_of_computing_square_roots
	private static final float invSqrt(float x)
	{
		float xhalf = 0.5f*x;
	        union { float x; int i; } u;
	
	        u.x = x;
	        u.i = 0x5f3759df - (u.i >> 1);
	        x = u.x * (1.5f - xhalf * u.x * u.x);
	        return x;
	}
	
	private static final float fastSqrt(float x)
	{
		return 1.0f / invSqrt(x);
	}
	*/
	
	private static final float sqrt(double value)
	{
		return (float)Math.sqrt(value);
	}

	private static final float cos(double angle)
	{
		return (float)Math.cos(angle);
	}

	private static final float sin(double angle)
	{
		return (float)Math.sin(angle);
	}
	
	private static final class SegmentData //TODO
	{
		int x0, y0, x1, y1, area, coll_count;
		
		public SegmentData(int x0, int y0, int x1, int y1, int area)
		{
			this.x0 = x0;
			this.y0 = y0;
			this.x1 = x1;
			this.y1 = y1;
			this.area = area;
			this.coll_count = 0;
		}
	}
	
	private final Random random;
	private int width, height, world_side;
	private float mass, velocity, alpha, vx, vy;

	float left;

	float top;

	private float cx, cy, dx, dy;
	private float[] map;
	private int[] age, segment;
	private List<SegmentData> seg_data = new ArrayList<SegmentData>();
	
	private int activeContinent;

	public final float getMomentum() { return mass * velocity; }
	public final int getHeight() { return height; }
	public final float getLeft() { return left; }
	public final float getTop() { return top; }
	public final float getVelocity() { return velocity; }
	public final float getVelX() { return vx; }
	public final float getVelY() { return vy; }
	public final int getWidth() { return width; }
	public final boolean isEmpty() { return mass <= 0; }
	
	public Plate(Random random, float[] m, int width, int height, int _x, int _y, int plate_age, int world_side)
	{
		this.random = random;
		this.width = width;
		this.height = height;
		this.world_side = world_side;
		this.mass = 0f;
		this.left = _x;
		this.top = _y;
		this.cx = this.cy = this.dx = this.dy = 0f;
		
		final int A = width * height; // A as in Area.
		final double angle = 2d * Math.PI * random.nextDouble();
		int i, j, k;
	
		if (m == null)
			return;
	
		map = new float[A];
		age = new int[A];
		segment = new int[A];
		
		for (int s = 0; s < segment.length; ++s)
			segment[s] = Integer.MAX_VALUE;
	
		velocity = 1;
		alpha = -(random.nextInt() & 1) * (float)Math.PI * 0.01f * random.nextFloat();
		vx = cos(angle) * INITIAL_SPEED_X;
		vy = sin(angle) * INITIAL_SPEED_X;
	
		for (j = k = 0; j < height; ++j)
			for (i = 0; i < width; ++i, ++k)
			{
				// Clone map data and count crust mass.
				mass += map[k] = m[k];
	
				// Calculate center coordinates weighted by mass.
				cx += i * m[k];
				cy += j * m[k];
	
				// Set the age of ALL points in this plate to same
				// value. The right thing to do would be to simulate
				// the generation of new oceanic crust as if the plate
				// had been moving to its current direction until all
				// plate's (oceanic) crust receive an age.
				age[k] = plate_age & -(m[k] > 0 ? 1 : 0);
			}
	
		// Normalize center of mass coordinates.
		cx /= mass;
		cy /= mass;
	}
	
	/*plate::~plate() throw()
	{
		delete[] map; map = 0;
		delete[] age; age = 0;
		delete[] segment; segment = 0;
	}*/
	
	private final void exit(int code)
	{
		System.exit(code);
	}
	
	public int addCollision(int wx, int wy)
	{
		MapIndex l = new MapIndex(wx, wy);
		int index = getMapIndex(l);
		int seg = seg_data.size();
	 
		if (DEBUG && index >= width * height)
		{
			System.out.println("Continental collision out of map bounds!");
			exit(1);
		}
	
		seg = segment[index];
	
		if (seg >= seg_data.size())
			seg = createSegment(l.x, l.y);
	
		if (DEBUG && seg >= seg_data.size())
		{
			System.out.println("Could not create segment!");
			exit(1);
		}
	
		++seg_data.get(seg).coll_count;
		return seg_data.get(seg).area;
	}
	
	public void addCrustByCollision(int x, int y, float z, int t)
	{
		// Add crust. Extend plate if necessary.
		setCrust(x, y, getCrust(x, y) + z, t);
	
		MapIndex m = new MapIndex(x, y);
		int index = getMapIndex(m);
		
		if (DEBUG && index >= width * height)
		{
			System.out.println("Aggregation went overboard!");
			exit(1);
		}
	
		segment[index] = activeContinent;
		SegmentData data = seg_data.get(activeContinent);
	
		++data.area;
		if (m.y < data.y0) data.y0 = m.y;
		if (m.y > data.y1) data.y1 = m.y;
		if (m.x < data.x0) data.x0 = m.x;
		if (m.x > data.x1) data.x1 = m.x;
	}


	/** Simulates subduction of oceanic plate under this plate.
	///
	/// Subduction is simulated by calculating the distance on surface
	/// that subducting sediment will travel under the plate until the
	/// subducting slab has reached certain depth where the heat triggers
	/// the melting and uprising of molten magma. 
	///
	/// @param	x	Origin of subduction on global world map (X).
	/// @param	y	Origin of subduction on global world map (Y).
	/// @param	z	Amount of sediment that subducts.
	/// @param	t	Time of creation of new crust.
	/// @param	dx	Direction of the subducting plate (X).
	/// @param	dy	Direction of the subducting plate (Y).*/
	public void addCrustBySubduction(int x, int y, float z, int t, float dx, float dy)
	{
		// TODO: Create an array of coordinate changes that would create
		//       a circle around current point. Array is static and it is
		//       initialized at the first call to this function.
		//       After all points of the circle are checked around subduction
		//       point the one with most land mass around it will be chosen as
		//       "most in land" point and subducting crust is added there.
		//       However to achieve a little more "natural" look normal
		//       distributed randomness is added around the "center" point.
		//       Benefits:
		//           NEVER adds crust aoutside plate.
		//           ALWAYS goes inland as much as possible
		//       Drawbacks:
		//           Additional logic required
		//           Might place crust on other continent on same plate!
		int index = getMapIndex(x, y);
	
		if (DEBUG && index >= width * height) // Should never be true!
		{
			System.out.println("Subduction origin not on plate!");
			System.out.println(String.format("%d, %d @ [%f, %f]x[%d, %d]\n", x, y, left, top, width, height));
			
			exit(1);
		}
	
		// Take vector difference only between plates that move more or less
		// to same direction. This makes subduction direction behave better.
		//
		// Use of "this" pointer is not necessary, but it make code clearer.
		// Cursed be those who use "m_" prefix in member names! >(
		
		float dot = this.vx * dx + this.vy * dy;
		dx -= this.vx * (dot > 0 ? 1 : 0);
		dy -= this.vy * (dot > 0 ? 1 : 0);
	
		float offset = random.nextFloat();
		offset *= offset * offset * (2 * (random.nextInt(2)) - 1);
		dx = 10 * dx + 3 * offset;
		dy = 10 * dx + 3 * offset;
	
		x = (int)((int)x + dx);
		y = (int)((int)y + dy);
	
		if (width == world_side)  x %= width;
		if (height == world_side) y %= height;
	
		index = y * width + x;
		if (index >= 0 && index < width * height && map[index] > 0) //TODO: remove neg check?
		{
			t = (age[index] + t) / 2;
			age[index] = t * (z > 0 ? 1 : 0);
	
			map[index] += z;
			mass += z;
		}
	}

	/** Add continental crust from this plate as part of other plate.
	///
	/// Aggregation of two continents is the event where the collided
	/// pieces of crust fuse together at the point of collision. It is
	/// crucial to merge not only the collided pieces of crust but also
	/// the entire continent that's part of the colliding tad of crust
	/// However, because one plate can contain many islands and pieces of
	/// continents, the merging must be done WITHOUT merging the entire
	/// plate and all those continental pieces that have NOTHING to do with
	/// the collision in question.
	///
	/// @param	p	Pointer to the receiving plate.
	/// @param	wx	X coordinate of collision point on world map.
	/// @param	wy	Y coordinate of collision point on world map.
	/// @return	Amount of crust aggregated to destination plate.
	 * */
	public float aggregateCrust(Plate p, int wx, int wy)
	{
		MapIndex l = new MapIndex(wx, wy);
		final int index = getMapIndex(l);
	
		if (DEBUG && index >= width * height)
		{
			System.out.println("Trying to aggregate beyond plate limits!");
			exit(1);
		}
	
		final int seg_id = segment[index];
	
		// This check forces the caller to do things in proper order!
		//
		// Usually continents collide at several locations simultaneously.
		// Thus if this segment that is being merged now is removed from
		// segmentation bookkeeping, then the next point of collision that is
		// processed during the same iteration step would cause the test
		// below to be true and system would experience a premature abort.
		//
		// Therefore, segmentation bookkeeping is left intact. It doesn't
		// cause significant problems because all crust is cleared and empty
		// points are not processed at all.

		if (DEBUG && seg_id >= seg_data.size())
		{
			System.out.println("Trying to aggregate without deforming first!");
			System.out.println(String.format("%d %d\n", wx, wy));
			exit(1);
		}
	
		// One continent may have many points of collision. If one of them
		// causes continent to aggregate then all successive collisions and
		// attempts of aggregation would necessarily change nothing at all,
		// because the continent was removed from this plate earlier!
		if (seg_data.get(seg_id).area == 0)
			return 0;	// Do not process empty continents.
	
		p.selectCollisionSegment(wx, wy);
	
		// Wrap coordinates around world edges to safeguard subtractions.
		wx += world_side;
		wy += world_side;
	
	//	printf("Aggregating segment [%d, %d]x[%d, %d] vs. [%d, %d]@[%d, %d]\n",
	//		seg_data[seg_id].x0, seg_data[seg_id].y0,
	//		seg_data[seg_id].x1, seg_data[seg_id].y1,
	//		width, height, lx, ly);
	
		float old_mass = mass;
	
		// Add all of the collided continent's crust to destination plate.
		for (int y = seg_data.get(seg_id).y0; y <= seg_data.get(seg_id).y1; ++y)
		  for (int x = seg_data.get(seg_id).x0; x <= seg_data.get(seg_id).x1; ++x)
		  {
			final int i = y * width + x;
			if ((segment[i] == seg_id) & (map[i] > 0))
			{
				p.addCrustByCollision(wx + x - l.x, wy + y - l.y,
					map[i], age[i]);
		
				mass -= map[i];
				map[i] = 0;
			}
		  }
	
		seg_data.get(seg_id).area = 0; // Mark segment as non-existent.
		return old_mass - mass;
	}
	
	public void applyFriction(float deformed_mass)
	{
		// Remove the energy that deformation consumed from plate's kinetic
		// energy: F - dF = ma - dF => a = dF/m.
		if (mass > 0)
		{
			float vel_dec = DEFORMATION_WEIGHT * deformed_mass / mass;
			vel_dec = vel_dec < velocity ? vel_dec : velocity;
	
			// Altering the source variable causes the order of calls to
			// this function to have difference when it shouldn't!
			// However, it's a hack well worth the outcome. :)
			velocity -= vel_dec;
		}
	}
	
	public void collide(Plate p, int wx, int wy, float coll_mass)
	{
		final float coeff_rest = 0.0f; // Coefficient of restitution.
		                              // 1 = fully elastic, 0 = stick together.
	
		// Calculate the normal to the curve/line at collision point.
		// The normal will point into plate B i.e. the "other" plate.
		//
		// Plates that wrap over world edges can mess the normal vector.
		// This could be solved by choosing the normal vector that points the
		// shortest path between mass centers but this causes problems when
		// plates are like heavy metal balls at a long rod and one plate's ball
		// collides at the further end of other plate's rod. Sure, this is
		// nearly never occurring situation but if we can easily do better then
		// why not do it?
		//
		// Better way is to select that normal vector that points along the
		// line that passes nearest the point of collision. Because point's
		// distance from line segment is relatively cumbersome to perform, the
		// vector is finally constructed as the sum of vectors <massCenterA, P> and
		// <P, massCenterB>. This solution works because collisions always
		// happen in the overlapping region of the two plates.
		MapIndex a = new MapIndex(wx, wy), b = new MapIndex(wx, wy);
		float ap_dx, ap_dy, bp_dx, bp_dy, nx, ny;
		int index = getMapIndex(a);
		int p_index = p.getMapIndex(b);
	
		if (DEBUG)
		if (index >= width * height || p_index >= p.width * p.height)
		{
			System.out.println(String.format("@%d, %d: out of colliding map's bounds!\n", wx, wy));
			exit(1);
		}

		ap_dx = (int)a.x - (int)cx;
		ap_dy = (int)a.y - (int)cy;
		bp_dx = (int)b.x - (int)p.cx;
		bp_dy = (int)b.y - (int)p.cy;
		nx = ap_dx - bp_dx;
		ny = ap_dy - bp_dy;
	
		if (nx * nx + ny * ny <= 0)
			return; // Avoid division by zero!
	
		// Scaling is required at last when impulses are added to plates!
		float n_len = sqrt(nx * nx + ny * ny);
		nx /= n_len;
		ny /= n_len;
	
		// Compute relative velocity between plates at the collision point.
		// Because torque is not included, calc simplifies to v_ab = v_a - v_b.
		final float rel_vx = vx - p.vx;
		final float rel_vy = vy - p.vy;
	
		// Get the dot product of relative velocity vector and collision vector.
		// Then get the projection of v_ab along collision vector.
		// Note that vector n must be a unit vector!
		final float rel_dot_n = rel_vx * nx + rel_vy * ny;
	
		if (rel_dot_n <= 0)
		{
	//		printf("n=%.2f, %.2f r=%.2f, %.2f, dot=%.4f\n",
	//			nx, ny, rel_vx, rel_vy, rel_dot_n);
			return; // Exit if objects are moving away from each other.
		}
	
		// Calculate the denominator of impulse: n . n * (1 / m_1 + 1 / m_2).
		// Use the mass of the colliding crust for the "donator" plate.
		float denom = (nx * nx + ny * ny) * (1.0f / mass + 1.0f / coll_mass);
	
		// Calculate force of impulse.
		float J = -(1 + coeff_rest) * rel_dot_n / denom;
	
		// Compute final change of trajectory.
		// The plate that is the "giver" of the impulse should receive a
		// force according to its pre-collision mass, not the current mass!
		dx += nx * J / mass;
		dy += ny * J / mass;
		p.dx -= nx * J / (coll_mass + p.mass);
		p.dy -= ny * J / (coll_mass + p.mass);
	
		// In order to prove that the code above works correctly, here is an
		// example calculation with ball A (mass 10) moving right at velocity
		// 1 and ball B (mass 100) moving up at velocity 1. Collision point
		// is at rightmost point of ball A and leftmost point of ball B.
		// Radius of both balls is 2.
		// ap_dx =  2;
		// ap_dy =  0;
		// bp_dx = -2;
		// bp_dy =  0;
		// nx = 2 - -2 = 4;
		// ny = 0 -  0 = 0;
		// n_len = sqrt(4 * 4 + 0) = 4;
		// nx = 4 / 4 = 1;
		// ny = 0 / 4 = 0;
		//
		// So far so good, right? Normal points into ball B like it should.
		//
		// rel_vx = 1 -  0 = 1;
		// rel_vy = 0 - -1 = 1;
		// rel_dot_n = 1 * 1 + 1 * 0 = 1;
		// denom = (1 * 1 + 0 * 0) * (1/10 + 1/100) = 1 * 11/100 = 11/100;
		// J = -(1 + 0) * 1 / (11/100) = -100/11;
		// dx = 1 * (-100/11) / 10 = -10/11;
		// dy = 0;
		// p.dx = -1 * (-100/11) / 100 = 1/11;
		// p.dy = -0;
		//
		// So finally:
		// vx = 1 - 10/11 = 1/11
		// vy = 0
		// p.vx = 0 + 1/11 = 1/11
		// p.vy = -1
		//
		// We see that in with restitution 0, both balls continue at same
		// speed along X axis. However at the same time ball B continues its
		// path upwards like it should. Seems correct right?
	}
	
	public void erode(float lower_bound)
	{
		float[] tmp = new float[width * height];
	
		mass = 0;
		cx = cy = 0;
	
		for (int y = 0; y < height; ++y)
	    for (int x = 0; x < width;  ++x)
	    {
			final int index = y * width + x;
			mass += map[index];
			tmp[index] += map[index]; // Careful not to overwrite earlier amounts.
		
			// Update the center coordinates weighted by mass.
			cx += x * map[index];
			cy += y * map[index];
		
			if (map[index] < lower_bound)
				continue;
		
			// Build masks for accessible directions (4-way).
			// Allow wrapping around map edges if plate has world wide dimensions.
			int w_mask = -((x > 0 ? 1 : 0) | (width == world_side ? 1 : 0));
			int e_mask = -((x < width - 1 ? 1 : 0) | (width == world_side ? 1 : 0));
			int n_mask = -((y > 0 ? 1 : 0) | (height == world_side ? 1 : 0));
			int s_mask = -((y < height - 1 ? 1 : 0) | (height == world_side ? 1 : 0));
		
			// Calculate the x and y offset of neighbor directions.
			// If neighbor is out of plate edges, set it to zero. This protects
			// map memory reads from segment faulting.
	    	int w = ((world_side + x - 1) % world_side) & w_mask;
	    	int e = ((world_side + x + 1) % world_side) & e_mask;
	    	int n = ((world_side + y - 1) % world_side) & n_mask;
	    	int s = ((world_side + y + 1) % world_side) & s_mask;
		
			// Calculate offsets within map memory.
			w = y * width + w;
			e = y * width + e;
			n = n * width + x;
			s = s * width + x;
		
			// Extract neighbors heights. Apply validity filtering: 0 is invalid.
			float w_crust = map[w] * (w_mask & (map[w] < map[index] ? 1 : 0));
			float e_crust = map[e] * (e_mask & (map[e] < map[index] ? 1 : 0));
			float n_crust = map[n] * (n_mask & (map[n] < map[index] ? 1 : 0));
			float s_crust = map[s] * (s_mask & (map[s] < map[index] ? 1 : 0));
		
			// Either this location has no neighbors (ARTIFACT!) or it is the
			// lowest part of its area. In either case the work here is done.
			if (w_crust + e_crust + n_crust + s_crust == 0)
				continue;
		
			// Calculate the difference in height between this point and its
			// nbours that are lower than this point.
			float w_diff = map[index] - w_crust;
			float e_diff = map[index] - e_crust;
			float n_diff = map[index] - n_crust;
			float s_diff = map[index] - s_crust;
		
			float min_diff = w_diff;
			min_diff -= (min_diff - e_diff) * (e_diff < min_diff ? 1 : 0);
			min_diff -= (min_diff - n_diff) * (n_diff < min_diff ? 1 : 0);
			min_diff -= (min_diff - s_diff) * (s_diff < min_diff ? 1 : 0);
		
			// Calculate the sum of difference between lower neighbors and
			// the TALLEST lower neighbor.
			float diff_sum = (w_diff - min_diff) * (w_crust > 0 ? 1 : 0) +
			                 (e_diff - min_diff) * (e_crust > 0 ? 1 : 0) +
			                 (n_diff - min_diff) * (n_crust > 0 ? 1 : 0) +
			                 (s_diff - min_diff) * (s_crust > 0 ? 1 : 0);
		
			if (DEBUG && diff_sum < 0)
			{
				System.out.println("Erosion difference sum is negative!");
				System.out.println(String.format("%f > %f %f %f %f\n", min_diff, w_diff, e_diff, n_diff, s_diff));
				exit(1);
			}
		
			if (diff_sum < min_diff)
			{
				// There's NOT enough room in neighbors to contain all the
				// crust from this peak so that it would be as tall as its
				// tallest lower neighbor. Thus first step is make ALL
				// lower neighbors and this point equally tall.
				tmp[w] += (w_diff - min_diff) * (w_crust > 0 ? 1 : 0);
				tmp[e] += (e_diff - min_diff) * (e_crust > 0 ? 1 : 0);
				tmp[n] += (n_diff - min_diff) * (n_crust > 0 ? 1 : 0);
				tmp[s] += (s_diff - min_diff) * (s_crust > 0 ? 1 : 0);
				tmp[index] -= min_diff;
		
				min_diff -= diff_sum;
		
				// Spread the remaining crust equally among all lower nbours.
				min_diff /= 1 
				+ (w_crust > 0 ? 1 : 0) 
				+ (e_crust > 0 ? 1 : 0) 
				+ (n_crust > 0 ? 1 : 0)
				+ (s_crust > 0 ? 1 : 0);
		
				tmp[w] += min_diff * (w_crust > 0 ? 1 : 0);
				tmp[e] += min_diff * (e_crust > 0 ? 1 : 0);
				tmp[n] += min_diff * (n_crust > 0 ? 1 : 0);
				tmp[s] += min_diff * (s_crust > 0 ? 1 : 0);
				tmp[index] += min_diff;
			}
			else
			{
				float unit = min_diff / diff_sum;
		
				// Remove all crust from this location making it as tall as
				// its tallest lower neighbour.
				tmp[index] -= min_diff;
		
				// Spread all removed crust among all other lower neighbours.
				tmp[w] += unit * (w_diff - min_diff) * (w_crust > 0 ? 1 : 0);
				tmp[e] += unit * (e_diff - min_diff) * (e_crust > 0 ? 1 : 0);
				tmp[n] += unit * (n_diff - min_diff) * (n_crust > 0 ? 1 : 0);
				tmp[s] += unit * (s_diff - min_diff) * (s_crust > 0 ? 1 : 0);
			}
	    }
	
	    //delete[] map;
	  	map = tmp;
	
		if (mass > 0)
		{
		    cx /= mass;
		    cy /= mass;
		}
	}
	
	public static final class CollisionInfo
	{
		public final int count;
		public final float ratio;
		
		private CollisionInfo(int count, float ratio)
		{
			this.count = count;
			this.ratio = ratio;
		}
	}
	
	public CollisionInfo getCollisionInfo(int wx, int wy) //TODO: make return value "lightweight"
	{
		int lx = wx, ly = wy;
		int index = getMapIndex(lx, ly);
		int seg = seg_data.size();
	
		if (DEBUG && index >= width * height)
		{
			System.out.println("getCollisionInfo: out of map bounds!");
			exit(1);
		}
	
		seg = segment[index];

		if (DEBUG && seg >= seg_data.size())
		{
			System.out.println("getCollisionInfo: no segment found!");
			exit(1);
		}
	
		return new CollisionInfo(
				seg_data.get(seg).coll_count,
				(float)seg_data.get(seg).coll_count / (float)(1 + seg_data.get(seg).area)); // +1 avoids DIV with zero.
	}
	
	public int getContinentArea(int wx, int wy)
	{
		final int index = getMapIndex(wx, wy);
	
		if (DEBUG && index >= width * height)
		{
			System.out.println("getContinentArea: out of map bounds!");
			exit(1);
		}
	
		if (DEBUG && segment[index] >= seg_data.size())
		{
			System.out.println("getContinentArea: no segment found!");
			exit(1);
		}
	
		return seg_data.get(segment[index]).area;
	}
	
	public float getCrust(int x, int y)
	{
		final int index = getMapIndex(x, y);
		return index < (int)(-1) ? map[index] : 0;
	}
	
	public int getCrustTimestamp(int x, int y)
	{
		final int index = getMapIndex(x, y);
		return index < (int)(-1) ? age[index] : 0;
	}
	
	/*public void getMap(float[] c, int[] t)
	{
		if (c) *c = map;
		if (t) *t = age;
	}*/
	
	public float[] getMap() { return map; }
	public int[] getAges() { return age; }
	
	@SuppressWarnings("unused")
	public void move()
	{
		float len;
	
		// Apply any new impulses to the plate's trajectory.
		vx += dx;
		vy += dy;
		dx = 0;
		dy = 0;
	
		// Force direction of plate to be unit vector.
		// Update velocity so that the distance of movement doesn't change.
		len = sqrt(vx*vx+vy*vy);
		vx /= len;
		vy /= len;
		velocity += len - 1.0;
		velocity *= velocity > 0 ? 1 : 0; // Round negative values to zero.
	
		// Apply some circular motion to the plate.
		if (false) // Or don't: it makes plates to halt sooner . more param tweaks.
		{
			float _cos = cos(alpha * velocity * velocity);
			float _sin = sin(alpha * velocity * velocity);
			float _vx = vx * _cos - vy * _sin;
			float _vy = vy * _cos + vx * _sin;
			vx = _vx;
			vy = _vy;
		}
	
		// Location modulations into range [0, world_side[ are a have to!
		// If left undone SOMETHING WILL BREAK DOWN SOMEWHERE in the code!
	
		if (DEBUG)
		if (left < 0 || left > world_side || top < 0 || top > world_side)
		{
			System.out.println("Location coordinates out of world map bounds (PRE)!");
			exit(1);
		}
	
		left += vx * velocity;
		left += left > 0 ? 0 : world_side;
		left -= left < world_side ? 0 : world_side;
	
		top += vy * velocity;
		top += top > 0 ? 0 : world_side;
		top -= top < world_side ? 0 : world_side;
	
		if (DEBUG)
		if (left < 0 || left > world_side || top < 0 || top > world_side)
		{
			System.out.println("Location coordinates out of world map bounds (POST)!");
			System.out.println(String.format("%f, %f, %f; %f, %f\n", vx, vy, velocity, left, top));
			exit(1);
		}
	}


	/** Clear any earlier continental crust partitions.
	///
	/// Plate has an internal bookkeeping of distinct areas of continental
	/// crust for more realistic collision response. However as the number
	/// of collisions that plate experiences grows, so does the bookkeeping
	/// of a continent become more and more inaccurate. Finally it results
	/// in striking artifacts that cannot overlooked.
	///
	/// To alleviate this problem without the need of per iteration
	/// recalculations plate supplies caller a method to reset its
	/// bookkeeping and start clean.*/
	public void resetSegments()
	{
		for (int i = 0; i < segment.length; ++i)
			segment[i] = Integer.MAX_VALUE;
		
		seg_data.clear();
	}

	/** Set the amount of plate's crustal material at some location.
	///
	/// If amount of crust to be set is negative, it'll be set to zero.
	///
	/// @param	x	Offset on the global world map along X axis.
	/// @param	y	Offset on the global world map along Y axis.
	/// @param	z	Amount of crust at given location.
	/// @param	t	Time of creation of new crust.
	 * */
	public void setCrust(int x, int y, float z, int t)
	{
		if (z < 0) // Do not accept negative values.
			z = 0;
	
		int index = getMapIndex(x, y);
	
		if (index >= width * height)
		{
			if (DEBUG && z <= 0)
			{
				System.out.println("Extending plate for nothing!");
				exit(1);
			}
	
			final int ilft = (int)left;
			final int itop = (int)top;
			final int irgt = ilft + width - 1;
			final int ibtm = itop + height - 1;
	
			x = (x + world_side) % world_side; // HACK!
			y = (y + world_side) % world_side; // Just to be safe...
	
			// Calculate distance of new point from plate edges.
			final int _lft = ilft - x;
			final int _rgt = (x < ilft ? world_side : 0) + x - irgt;
			final int _top = itop - y;
			final int _btm = (y < itop ? world_side : 0) + y - ibtm;
	
			// Set larger of horizontal/vertical distance to zero.
			// A valid distance is NEVER larger than world's side's length!
			int d_lft = _lft & -(_lft <  _rgt ? 1 : 0) & -(_lft < world_side ? 1 : 0);
			int d_rgt = _rgt & -(_rgt <= _lft ? 1 : 0) & -(_rgt < world_side ? 1 : 0);
			int d_top = _top & -(_top <  _btm ? 1 : 0) & -(_top < world_side ? 1 : 0);
			int d_btm = _btm & -(_btm <= _top ? 1 : 0) & -(_btm < world_side ? 1 : 0);
	
			// Scale all changes to multiple of 8.
			d_lft = ((d_lft > 0 ? 1 : 0) + (d_lft >> 3)) << 3;
			d_rgt = ((d_rgt > 0 ? 1 : 0) + (d_rgt >> 3)) << 3;
			d_top = ((d_top > 0 ? 1 : 0) + (d_top >> 3)) << 3;
			d_btm = ((d_btm > 0 ? 1 : 0) + (d_btm >> 3)) << 3;
	
			// Make sure plate doesn't grow bigger than the system it's in!
			if (width + d_lft + d_rgt > world_side)
			{
				d_lft = 0;
				d_rgt = world_side - width;
			}
	
			if (height + d_top + d_btm > world_side)
			{
				d_top = 0;
				d_btm = world_side - height;
			}
	
			if (DEBUG)
			if (d_lft + d_rgt + d_top + d_btm == 0)
			{
				System.out.println(String.format("[%d, %d]x[%d, %d], [%d, %d]/[%d, %d]\n",
					(int)left, (int)top, (int)left+width,
					(int)top+height,
					x + world_side * (x < world_side ? 1 : 0),
					y + world_side * (y < world_side ? 1 : 0),
					x % world_side, y % world_side));
	
				System.out.println("Index out of bounds, but nowhere to grow!");
				exit(1);
			}
	
			final int old_width = width;
			final int old_height = height;
			
			left -= d_lft;
			left += left >= 0 ? 0 : world_side;
			width += d_lft + d_rgt;
	
			top -= d_top;
			top += top >= 0 ? 0 : world_side;
			height += d_top + d_btm;
	
	//		printf("%dx%d + [%d,%d] + [%d, %d] = %dx%d\n",
	//			old_width, old_height,
	//			d_lft, d_top, d_rgt, d_btm, width, height);
	
			float[] tmph = new float[width * height];
			int[] tmpa = new int[width * height];
			int[] tmps = new int[width * height];
			
			for (int i = 0; i < tmps.length; ++i)
				tmps[i] = Integer.MAX_VALUE;
	
			// copy old plate into new.
			for (int j = 0; j < old_height; ++j)
			{
				final int dest_i = (d_top + j) * width + d_lft;
				final int src_i = j * old_width;
				
				for (int c = 0; c < old_width; ++c)
				{
					tmph[dest_i + c] = map[src_i + c];
					tmpa[dest_i + c] = age[src_i + c];
					tmps[dest_i + c] = segment[src_i + c];
				}
			}
	
			//delete[] map;
			//delete[] age;
			//delete[] segment;
			map = tmph;
			age = tmpa;
			segment = tmps;
	
			// Shift all segment data to match new coordinates.
			for (int s = 0; s < seg_data.size(); ++s)
			{
				seg_data.get(s).x0 += d_lft;
				seg_data.get(s).x1 += d_lft;
				seg_data.get(s).y0 += d_top;
				seg_data.get(s).y1 += d_top;
			}
	
			index = getMapIndex(x, y);
	
			if (DEBUG && index >= width * height)
			{
				System.out.println(String.format("Index out of bounds after resize!\n" + 
					"[%d, %d]x[%d, %d], [%d, %d]/[%d, %d]\n",
					(int)left, (int)top, (int)left+width,
					(int)top+height,
					x, y, x % world_side, y % world_side));
				exit(1);
			}
		}
	
		// Update crust's age.
		// If old crust exists, new age is mean of original and supplied ages.
		// If no new crust is added, original time remains intact.
		final int old_crust = -(map[index] > 0 ? 1 : 0);
		final int new_crust = -(z > 0 ? 1 : 0);
		t = (t & ~old_crust) | ((age[index] + t) / 2 & old_crust);
		age[index] = (t & new_crust) | (age[index] & ~new_crust);
	
		mass -= map[index];
		map[index] = z;		// Set new crust height to desired location.
		mass += z;		// Update mass counter.
	}

	/** Remember the currently processed continent's segment number.
	///
	/// @param	coll_x	Origin of collision on global world map (X).
	/// @param	coll_y	Origin of collision on global world map (Y).
	*/
	public void selectCollisionSegment(int coll_x, int coll_y)
	{
		int index = getMapIndex(coll_x, coll_y);
	
		activeContinent = seg_data.size();

		if (DEBUG && index >= width * height)
		{
			System.out.println("Collision segment cannot be set outside plate!");
			exit(1);
		}
	
		activeContinent = segment[index];
	
		if (DEBUG && activeContinent >= seg_data.size())
		{
			System.out.println("Collision happened at unsegmented location!");
			exit(1);
		}
	}
	
	///////////////////////////////////////////////////////////////////////////////
	/// Private methods ///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////

	/** Separate a continent at (X, Y) to its own partition.
	 * 
	 * Method analyzes the pixels 4-ways adjacent at the given location
	 * and labels all connected continental points with same segment ID.
	 * 
	 * @param	x	Offset on the local height map along X axis.
	 * @param	y	Offset on the local height map along Y axis.
	 *
	 * @return	ID of created segment on success, otherwise -1.
	 * */
	private int createSegment(int x, int y)
	{
		final int origin_index = y * width + x;
		final int ID = seg_data.size();
	
		if (segment[origin_index] < ID)
			return segment[origin_index];
	
		boolean canGoLeft =  x > 0          && map[origin_index - 1] >= CONT_BASE;
		boolean canGoRight = x < width - 1  && map[origin_index+1] >= CONT_BASE;
		boolean canGoUp =    y > 0          && map[origin_index - width] >= CONT_BASE;
		boolean canGoDown =  y < height - 1 && map[origin_index + width] >= CONT_BASE;
		int nbour_id = ID;
	
		// This point belongs to no segment yet.
		// However it might be a neighbour to some segment created earlier.
		// If such neighbour is found, associate this point with it.
		if (canGoLeft && segment[origin_index - 1] < ID)
			nbour_id = segment[origin_index - 1];
		else if (canGoRight && segment[origin_index + 1] < ID)
			nbour_id = segment[origin_index + 1];
		else if (canGoUp && segment[origin_index - width] < ID)
			nbour_id = segment[origin_index - width];
		else if (canGoDown && segment[origin_index + width] < ID)
			nbour_id = segment[origin_index + width];
	
		if (nbour_id < ID)
		{
			segment[origin_index] = nbour_id;
			seg_data.get(nbour_id).area += 1;
	
			if (y < seg_data.get(nbour_id).y0) seg_data.get(nbour_id).y0 = y;
			if (y > seg_data.get(nbour_id).y1) seg_data.get(nbour_id).y1 = y;
			if (x < seg_data.get(nbour_id).x0) seg_data.get(nbour_id).x0 = x;
			if (x > seg_data.get(nbour_id).x1) seg_data.get(nbour_id).x1 = x;
	
			return nbour_id;
		}
	
		int lines_processed;
		SegmentData data = new SegmentData(x, y, x, y, 0);
	
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] spans_todo = new ArrayList[height];
		@SuppressWarnings("unchecked")
		ArrayList<Integer>[] spans_done = new ArrayList[height];
		
		for (int h = 0; h < height; ++h)
		{
			spans_todo[h] = new ArrayList<Integer>();
			spans_done[h] = new ArrayList<Integer>();
		}
	
		segment[origin_index] = ID;
		spans_todo[y].add(x);
		spans_todo[y].add(x);
	
		do
		{
			lines_processed = 0;
			for (int line = 0; line < height; ++line)
			{
				int start, end;
		
				if (spans_todo[line].size() == 0)
					continue;
		
				do // Find an unscanned span on this line.
				{
					end = spans_todo[line].get(spans_todo[line].size() - 1);
					spans_todo[line].remove(spans_todo[line].size() - 1);
	
					start = spans_todo[line].get(spans_todo[line].size() - 1);
					spans_todo[line].remove(spans_todo[line].size() - 1);
		
					// Reduce any done spans from this span.
					for (int j = 0; j < spans_done[line].size(); j += 2)
					{
						// Saved coordinates are AT the point
						// that was included last to the span.
						// That's why equalities matter.
		
						if (start >= spans_done[line].get(j) &&
						    start <= spans_done[line].get(j+1))
							start = spans_done[line].get(j+1) + 1;
		
						if (end >= spans_done[line].get(j) &&
						    end <= spans_done[line].get(j+1))
							end = spans_done[line].get(j) - 1;
					}
		
					// Unsigned-ness hacking!
					// Required to fix the underflow of end - 1.
					start |= -(end >= width ? 1 : 0);
					end -= (end >= width ? 1 : 0);
		
				}
				while (start > end && !spans_todo[line].isEmpty());
		
				if (start > end) // Nothing to do here anymore...
					continue;
	
				// Calculate line indices. Allow wrapping around map edges.
				final int row_above = ((line - 1) & -(line > 0 ? 1 : 0)) |
					((height - 1) & -(line == 0 ? 1 : 0));
				final int row_below = (line + 1) & -(line < height - 1 ? 1 : 0);
				final int line_here = line * width;
				final int line_above = row_above * width;
				final int line_below = row_below * width;
		
				// Extend the beginning of line.
				while (start > 0 && segment[line_here+start-1] > ID &&
					map[line_here+start-1] >= CONT_BASE)
				{
					--start;
					segment[line_here + start] = ID;
		
					// Count volume of pixel...
				}
		
				// Extend the end of line.
				while (end < width - 1 &&
					segment[line_here + end + 1] > ID &&
					map[line_here + end + 1] >= CONT_BASE)
				{
					++end;
					segment[line_here + end] = ID;
		
					// Count volume of pixel...
				}
		
				// Check if should wrap around left edge.
				if (width == world_side && start == 0 &&
					segment[line_here+width-1] > ID &&
					map[line_here+width-1] >= CONT_BASE)
				{
					segment[line_here + width - 1] = ID;
					spans_todo[line].add(width - 1);
					spans_todo[line].add(width - 1);
		
					// Count volume of pixel...
				}
		
				// Check if should wrap around right edge.
				if (width == world_side && end == width - 1 &&
					segment[line_here+0] > ID &&
					map[line_here+0] >= CONT_BASE)
				{
					segment[line_here + 0] = ID;
					spans_todo[line].add(0);
					spans_todo[line].add(0);
		
					// Count volume of pixel...
				}
		
				data.area += 1 + end - start; // Update segment area counter.
		
				// Record any changes in extreme dimensions.
				if (line < data.y0) data.y0 = line;
				if (line > data.y1) data.y1 = line;
				if (start < data.x0) data.x0 = start;
				if (end > data.x1) data.x1 = end;
		
				if (line > 0 || height == world_side)
				for (int j = start; j <= end; ++j)
					if (segment[line_above + j] > ID && map[line_above + j] >= CONT_BASE)
					{
						int a = j;
						segment[line_above + a] = ID;
			
						// Count volume of pixel...
			
						while (++j < width &&
						       segment[line_above + j] > ID &&
						       map[line_above + j] >= CONT_BASE)
						{
							segment[line_above + j] = ID;
			
							// Count volume of pixel...
						}
			
						int b = --j; // Last point is invalid.
			
						spans_todo[row_above].add(a);
						spans_todo[row_above].add(b);
						++j; // Skip the last scanned point.
					}
	
				if (line < height - 1 || height == world_side)
				for (int j = start; j <= end; ++j)
				    if (segment[line_below + j] > ID &&  map[line_below + j] >= CONT_BASE)
				    {
						int a = j;
						segment[line_below + a] = ID;
			
						// Count volume of pixel...
			
						while (++j < width &&
						       segment[line_below + j] > ID &&
						       map[line_below + j] >= CONT_BASE)
						{
							segment[line_below + j] = ID;
			
							// Count volume of pixel...
						}
			
						int b = --j; // Last point is invalid.
			
						spans_todo[row_below].add(a);
						spans_todo[row_below].add(b);
						++j; // Skip the last scanned point.
				    }
	
				spans_done[line].add(start);
				spans_done[line].add(end);
				++lines_processed;
			}
		}
		while (lines_processed > 0);
	
		seg_data.add(data);
	//	printf("Created segment [%d, %d]x[%d, %d]@[%d, %d].\n",
	//		data.x0, data.y0, data.x1, data.y1, x, y);
	
		return ID;
	}

	public static final class MapIndex
	{
		public int x, y;
		
		public MapIndex(int x, int y)
		{
			this.x = x;
			this.y = y;
		}
	}

	private final MapIndex temp = new MapIndex(0, 0);
	
	public int getMapIndex(int x, int y)
	{
		temp.x = x;
		temp.y = y;
		return getMapIndex(temp);
	}

	/** 
	 * Get pointers to plate's data.
	 * @param	c	Address of crust height map is stored here.
	 * @param	t	Address of crust timestamp map is stored here.
	 * */
	public int getMapIndex(MapIndex mapIndex)
	{
		int x = mapIndex.x;
		int y = mapIndex.y;
		final int ilft = (int)left;
		final int itop = (int)top;
		final int irgt = ilft + width;
		final int ibtm = itop + height;

		//TODO: "+ world_side" might be unnecessary
		x = (x + world_side) % world_side; // Sometimes input is beyond map dimensions.
		y = (y + world_side) % world_side; // Scale it to fit within world map.

		///////////////////////////////////////////////////////////////////////
		// If you think you're smart enough to optimize this then PREPARE to be
		// smart as HELL to debug it!
		///////////////////////////////////////////////////////////////////////

		final boolean xOkA = x >= ilft && x < irgt;
		final boolean xOkB = x + world_side >= ilft && x + world_side < irgt;
		final boolean xOk = xOkA || xOkB;

		final boolean yOkA = y >= itop && y < ibtm;
		final boolean yOkB = y + world_side >= itop && y + world_side < ibtm;
		final boolean yOk = yOkA || yOkB;

		x += world_side & -(x < ilft ? 1 : 0); // Point is within plate's map: wrap
		y += world_side & -(y < itop ? 1 : 0); // it around world edges if necessary.

		x -= ilft; // Calculate offset within local map.
		y -= itop;

		int failMask = (!xOk || !yOk) ? -1 : 0;

		if (DEBUG)
		{
			if (failMask != 0)
			{
				boolean X_OK = (mapIndex.x >= ilft && mapIndex.x < irgt) || 
					(mapIndex.x + world_side >= ilft && mapIndex.x + world_side < irgt);
				boolean Y_OK = (mapIndex.y >= itop && mapIndex.y < ibtm) || 
					(mapIndex.y + world_side >= itop && mapIndex.y + world_side < ibtm);
	
				if (X_OK && Y_OK)
				{
					System.out.println("MapIndex has an error, goddamn!");
					exit(1);
				}
			}
			else if (y >= height || x >= width)
			{
				System.out.println(String.format("\nMap Index error:\n%d <= %d < %d, %d <= %d < %d",
					0, x, width, 0, y, height));
				System.out.println(String.format("%d <= %d < %d, %d <= %d < %d",
						ilft, mapIndex.x, irgt, itop, mapIndex.y, ibtm));
				exit(1);
			}
		}

		mapIndex.x = x & ~failMask;
		mapIndex.y = y & ~failMask;
		
		return (y * width + x) | (failMask & Integer.MAX_VALUE);
	}
}

