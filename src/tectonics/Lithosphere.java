package tectonics;

import java.util.ArrayList;
import java.util.Random;

import src.Array2D;
import src.IModifier;
import tectonics.Plate.CollisionInfo;

public class Lithosphere implements IModifier<Array2D>
{
	/**
	 * Wrapper for growing plate from a seed. Contains plate's dimensions.
	 *
	 * Used exclusively in plate creation phase.
	 */
	public static final class PlateArea
	{
		public ArrayList<Integer> border; ///< Plate's unprocessed border pixels.
		public int btm; ///< Most bottom pixel of plate.
		public int lft; ///< Most left pixel of plate.
		public int rgt; ///< Most right pixel of plate.
		public int top; ///< Most top pixel of plate.
		public int wdt; ///< Width of area in pixels.
		public int hgt; ///< Height of area in pixels.
	}


	/**
	 * Container for collision details between two plates.
	 *
	 * In simulation there's usually 2-5 % collisions of the entire map
	 * area. In a 512*512 map that means 5000-13000 collisions.
	 *
	 * When plate collisions are recorded and processed pair-by-pair, some
	 * of the information is lost if more than two plates collide at the
	 * same point (there will be no record of the two lower plates
	 * colliding together, just that they both collided with the tallest
	 * plate) ONLY IF ALL the collisions between ANY TWO plates of that
	 * group always include a third, taller/higher  plate. This happens
	 * most often when plates have long, sharp spikes i.e. in the
	 * beginning.
	 */
	class PlateCollision
	{
		public PlateCollision(int index, int wx, int wy, float crust)
		{
			this.index = index;
			this.wx = wx;
			this.wy = wy;
			this.crust = crust;
		}

		public final int index; ///< Index of the other plate involved in the event.
		public final int wx, wy; ///< Coordinates of collision in world space.
		public final float crust; ///< Amount of crust that will deform/subduct.
	}
	
	private static final void exit(int code)
	{
		System.exit(code);
	}
	
	private static final void puts(String str, Object... toFormat)
	{
		System.out.println(String.format(str, toFormat));
	}
	
	@SuppressWarnings("unused")
	private static final boolean isPowerOfTwo(int n)
	{
		return ((n != 0) && (n & (n - 1)) == 0);
	}
	
	private static boolean DEBUG = true;

	private static final boolean BOOL_REGENERATE_CRUST = true;
	
	private static final float BUOYANCY_BONUS_X = 3;
	private static final int MAX_BUOYANCY_AGE = 20;
	private static final float MULINV_MAX_BUOYANCY_AGE = 1.0f / (float)MAX_BUOYANCY_AGE;
	
	private static final float RESTART_ENERGY_RATIO = 0.15f;
	private static final float RESTART_SPEED_LIMIT = 2.0f;
	private static final int NO_COLLISION_TIME_LIMIT = 10;

	private static final float CONTINENTAL_BASE = 1.0f;
	private static final float OCEANIC_BASE     = 0.1f;
	
	private static final float FLOAT_EPSILON = 1E-8f;
	
	private final Random random;
	private final float sea_level, folding_ratio, aggr_overlap_rel;
	private final int num_plates, aggr_overlap_abs, erosion_period, max_cycles;
	
	private Array2D hmap;
	private int[] imap;
	private Plate[] plates;
	private int map_side, iter_count, last_coll_count, cycle_count;
	private float peak_Ek;
	
	private ArrayList<ArrayList<PlateCollision>> collisions = new ArrayList<ArrayList<PlateCollision>>();
	private ArrayList<ArrayList<PlateCollision>> subductions = new ArrayList<ArrayList<PlateCollision>>();
	
	//int findBound(final int[] map, int length, int x0, int y0, int dx, int dy);
	
	//int findPlate(Plate[][] plates, float x, float y, int num_plates);
	
	public Lithosphere(Random random, int num_plates, float sea_level,
		int erosion_period, float folding_ratio, int aggr_ratio_abs,
		float aggr_ratio_rel, int max_cycles)
		/*hmap(0), plates(0), aggr_overlap_abs(aggr_ratio_abs),
		aggr_overlap_rel(aggr_ratio_rel), cycle_count(0),
		erosion_period(_erosion_period), folding_ratio(_folding_ratio),
		iter_count(0), map_side(map_side_length), max_cycles(num_cycles),
		num_plates(0)*/
	{
		this.random = random;
		this.num_plates = num_plates;
		this.sea_level = sea_level;
		this.erosion_period = erosion_period;
		this.folding_ratio = folding_ratio;
		this.aggr_overlap_abs = aggr_ratio_abs;
		this.aggr_overlap_rel = aggr_ratio_rel;
		this.max_cycles = max_cycles;
	}
	
	/*~lithosphere()
	{
		delete[] plates; plates = 0;
		delete[] imap;   imap = 0;
		delete[] hmap;   hmap = 0;
	}*/
	
	@Override
	public void modify(Array2D tmp)
	{
		map_side = tmp.width;
		final int A = map_side * map_side;
		
		if (tmp.width != tmp.height)
			throw new Error("heightMap width and height must be equal!");
	
		float lowest = tmp.get1D(0), highest = tmp.get1D(0);
		for (int i = 1; i < A; ++i)
		{
			lowest = lowest < tmp.get1D(i) ? lowest : tmp.get1D(i);
			highest = highest > tmp.get1D(i) ? highest : tmp.get1D(i);
		}
	
		for (int i = 0; i < A; ++i) // Scale to [0 ... 1]
			tmp.set1D(i, (tmp.get1D(i) - lowest) / (highest - lowest));
	
		float sea_threshold = 0.5f;
		float th_step = 0.5f;
	
		// Find the actual value in height map that produces the continent-sea
		// ratio defined be "sea_level".
		while (th_step > 0.01)
		{
			int count = 0;
			for (int i = 0; i < A; ++i)
				count += (tmp.get1D(i) < sea_threshold ? 1 : 0);
	
			th_step *= 0.5;
			if (count / (float)A < sea_level)
				sea_threshold += th_step;
			else
				sea_threshold -= th_step;
		}
	
		for (int i = 0; i < A; ++i) // Genesis 1:9-10.
		{
	//		if (tmp[i] < sea_level)
	//			tmp[i] /= sea_level;
	//		else
	//			tmp[i] = 1.0+(tmp[i] - sea_level) * 20;
			tmp.set1D(i, (tmp.get1D(i) > sea_threshold ? 1 : 0) *
				(tmp.get1D(i) + CONTINENTAL_BASE) + 
				(tmp.get1D(i) <= sea_threshold ? 1 : 0) * OCEANIC_BASE);
		}
		
		hmap = tmp;
		imap = new int[map_side * map_side];
		
		createPlates();
		
		for (int i = 0; i < 600 * 2; ++i)
			update();
	}
	
	private void createPlates()
	{
		final int map_area = map_side * map_side;
	
		ArrayList<PlateCollision> vec = new ArrayList<PlateCollision>();
		vec.ensureCapacity(map_side * 4); // == map's circumference.
	
		collisions.ensureCapacity(collisions.size() + num_plates);
		subductions.ensureCapacity(subductions.size() + num_plates);
	
		for (int i = 0; i < num_plates; ++i)
		{
			collisions.add(vec);
			subductions.add(vec);
		}
	
		// Initialize "Free plate center position" lookup table.
		// This way two plate centers will never be identical.
		for (int i = 0; i < map_area; ++i)
			imap[i] = i;
	
		// Select N plate centers from the global map.
		PlateArea[] area = new PlateArea[num_plates];
		for (int i = 0; i < num_plates; ++i)
		{
			// Randomly select an unused plate origin.
			final int p = imap[random.nextInt(map_area - i)];
			final int y = p / map_side;
			final int x = p - y * map_side;

			area[i] = new PlateArea();
			
			area[i].lft = area[i].rgt = x; // Save origin...
			area[i].top = area[i].btm = y;
			area[i].wdt = area[i].hgt = 1;
	
			area[i].border = new ArrayList<Integer>(8);
			area[i].border.add(p); // ...and mark it as border.
	
			// Overwrite used entry with last unused entry in array.
			imap[p] = imap[map_area - i - 1];
		}
	
		final int[] owner = imap; // Create an alias.
		for (int i = 0; i < map_area; ++i)
			owner[i] = Integer.MAX_VALUE;
	
		// "Grow" plates from their origins until surface is fully populated.
		int max_border = 1;
		int i;
		while (max_border != 0)
			for (max_border = i = 0; i < num_plates; ++i)
			{
				final int N = area[i].border.size();
				max_border = max_border > N ? max_border : N;
	
				if (N == 0)
					continue;
	
				final int j = random.nextInt(N);
				final int p = area[i].border.get(j);
				final int cy = p / map_side;
				final int cx = p - cy * map_side;
	
				final int lft = cx > 0 ? cx - 1 : map_side - 1;
				final int rgt = cx < map_side - 1 ? cx + 1 : 0;
				final int top = cy > 0 ? cy - 1 : map_side - 1;
				final int btm = cy < map_side - 1 ? cy + 1 : 0;
	
				final int n = top * map_side +  cx; // North.
				final int s = btm * map_side +  cx; // South.
				final int w =  cy * map_side + lft; // West.
				final int e =  cy * map_side + rgt; // East.
	
				if (owner[n] >= num_plates)
				{
					owner[n] = i;
					area[i].border.add(n);
	
					if (area[i].top == ((top + 1) % map_side))
					{
						area[i].top = top;
						area[i].hgt++;
					}
				}
	
				if (owner[s] >= num_plates)
				{
					owner[s] = i;
					area[i].border.add(s);
	
					if (btm == ((area[i].btm + 1) % map_side))
					{
						area[i].btm = btm;
						area[i].hgt++;
					}
				}
	
				if (owner[w] >= num_plates)
				{
					owner[w] = i;
					area[i].border.add(w);
	
					if (area[i].lft == ((lft + 1) % map_side))
					{
						area[i].lft = lft;
						area[i].wdt++;
					}
				}
	
				if (owner[e] >= num_plates)
				{
					owner[e] = i;
					area[i].border.add(e);
	
					if (rgt == ((area[i].rgt + 1) % map_side))
					{
						area[i].rgt = rgt;
						area[i].wdt++;
					}
				}
	
				// Overwrite processed point with unprocessed one.
				area[i].border.set(j, area[i].border.get(area[i].border.size() - 1));
				area[i].border.remove(area[i].border.size() - 1);
			}
	
		plates = new Plate[num_plates];
	
		// Extract and create plates from initial terrain.
		for (i = 0; i < num_plates; ++i)
		{
			area[i].wdt = area[i].wdt < map_side ? area[i].wdt : map_side - 1;
			area[i].hgt = area[i].hgt < map_side ? area[i].hgt : map_side - 1;
	
			final int x0 = area[i].lft;
			final int x1 = 1 + x0 + area[i].wdt;
			final int y0 = area[i].top;
			final int y1 = 1 + y0 + area[i].hgt;
			final int width = x1 - x0;
			final int height = y1 - y0;
			float[] plt = new float[width * height];
	
	//		printf("plate %d: (%d, %d)x(%d, %d)\n", i, x0, y0, width,
	//			height);
	
			// Copy plate's height data from global map into local map.
			for (int y = y0, j = 0; y < y1; ++y)
				for (int x = x0; x < x1; ++x, ++j)
				{
					int k = (y % map_side) * map_side + (x % map_side);
					plt[j] = hmap.get1D(k) * (owner[k] == i ? 1 : 0);
				}
	
			// Create plate.
			plates[i] = new Plate(random, plt, width, height, x0, y0, i, map_side);
			//delete[] plt;
		}
	
		iter_count = num_plates + MAX_BUOYANCY_AGE;
		peak_Ek = 0;
		last_coll_count = 0;
		//delete[] area;
	}
	
	public final int getPlateCount()
	{
		return num_plates;
	}
	
	@SuppressWarnings("unused")
	public void update()
	{
		float totalVelocity = 0;
		float systemKineticEnergy = 0;
	
		for (int i = 0; i < num_plates; ++i)
		{
			totalVelocity += plates[i].getVelocity();
			systemKineticEnergy += plates[i].getMomentum();
		}
	
		if (systemKineticEnergy > peak_Ek)
			peak_Ek = systemKineticEnergy;
	
	//	printf("%f > %f, ", totalVelocity, RESTART_SPEED_LIMIT);
	//	printf("%f/%f = %f > %f\n", systemKineticEnergy, peak_Ek,
	//		systemKineticEnergy / peak_Ek, RESTART_ENERGY_RATIO);
	
		// If there's no continental collisions during past iterations,
		// then interesting activity has ceased and we should restart.
		// Also if the simulation has been going on for too long already,
		// restart, because interesting stuff has most likely ended.
		if (totalVelocity < RESTART_SPEED_LIMIT ||
		    systemKineticEnergy / peak_Ek < RESTART_ENERGY_RATIO ||
		    last_coll_count > NO_COLLISION_TIME_LIMIT ||
		    iter_count > 600)
		{
			restart();
			return;
		}
	
		final int map_area = map_side * map_side;
		final int[] prev_imap = imap;
		int[] amap = new int[map_area];
		imap = new int[map_area];
	
		// Realize accumulated external forces to each plate.
		for (int i = 0; i < num_plates; ++i)
		{
			// Dont't do it yet.. There's problems with the index map...
			/*if (plates[i].isEmpty())
			{
				//delete plates[i];
				plates[i] = plates[--num_plates];
				--i;
	
				continue;
			}*/
	
			plates[i].resetSegments();
	
			if (erosion_period > 0 && iter_count % erosion_period == 0)
				plates[i].erode(CONTINENTAL_BASE);
	
			plates[i].move();
		}
	
	//	static int max_collisions = 0;	// DEBUG!!!
		int oceanic_collisions = 0;
		int continental_collisions = 0;
	
		// Update height and plate index maps.
		// Doing it plate by plate is much faster than doing it index wise:
		// Each plate's map's memory area is accessed sequentially and only
		// once as opposed to calculating "num_plates" indices within plate
		// maps in order to find out which plate(s) own current location.
		
		for (int i = 0; i < map_area; ++i)
		{
			hmap.set1D(i, 0f);
			imap[i] = Integer.MAX_VALUE;
		}
		
		for (int i = 0; i < num_plates; ++i)
		{
			final int x0 = (int)plates[i].getLeft();
			final int y0 = (int)plates[i].getTop();
			final int x1 = x0 + plates[i].getWidth();
			final int y1 = y0 + plates[i].getHeight();
		
			final float[] this_map = plates[i].getMap();
			final int[] this_age = plates[i].getAges();
		
			// Copy first part of plate onto world map.
			for (int y = y0, j = 0; y < y1; ++y)
		    for (int x = x0;        x < x1; ++x, ++j)
		    {
				final int x_mod = x % map_side;
				final int y_mod = y % map_side;
		
				final int k = y_mod * map_side + x_mod;
		
				if (this_map[j] < 2 * FLOAT_EPSILON) // No crust here...
					continue;
		
				if (imap[k] >= num_plates) // No one here yet?
				{
					// This plate becomes the "owner" of current location
					// if it is the first plate to have crust on it.
					hmap.set1D(k, this_map[j]);
					imap[k] = i;
					amap[k] = this_age[j];
		
					continue;
				}
	
				// DO NOT ACCEPT HEIGHT EQUALITY! Equality leads to subduction
				// of shore that 's barely above sea level. It's a lot less
				// serious problem to treat very shallow waters as continent...
				final boolean prev_is_oceanic = hmap.get1D(k) < CONTINENTAL_BASE;
				final boolean this_is_oceanic = this_map[j] < CONTINENTAL_BASE;
		
				final int prev_timestamp = plates[imap[k]].getCrustTimestamp(x_mod, y_mod);
				final int this_timestamp = this_age[j];
				
				final boolean prev_is_bouyant = hmap.get1D(k) > this_map[j] ||
					(hmap.get1D(k) + 2 * FLOAT_EPSILON > this_map[j] &&
					 hmap.get1D(k) < 2 * FLOAT_EPSILON + this_map[j] &&
					 prev_timestamp >= this_timestamp);
	
				// Handle subduction of oceanic crust as special case.
				if (this_is_oceanic && prev_is_bouyant)
				{
					// This plate will be the subducting one.
					// The level of effect that subduction has
					// is directly related to the amount of water
					// on top of the subducting plate.
					final float sediment = OCEANIC_BASE * (1f - this_map[j] / CONTINENTAL_BASE);
		
					// Save collision to the receiving plate's list.
					PlateCollision coll = new PlateCollision(i, x_mod, y_mod, sediment);
					subductions.get(imap[k]).add(coll);
					++oceanic_collisions;
		
					// Remove subducted oceanic lithosphere from plate.
					// This is crucial for
					// a) having correct amount of colliding crust (below)
					// b) protecting subducted locations from receiving
					//    crust from other subductions/collisions.
					plates[i].setCrust(
							x_mod,
							y_mod, 
							this_map[j] - OCEANIC_BASE,
							this_timestamp);
		
					if (this_map[j] <= 0)
						continue; // Nothing more to collide.
				}
				else if (prev_is_oceanic)
				{
					final float sediment = OCEANIC_BASE *
						(CONTINENTAL_BASE - hmap.get1D(k)) /
						CONTINENTAL_BASE;
		
					PlateCollision coll = new PlateCollision(imap[k], x_mod, y_mod, sediment);
					subductions.get(i).add(coll);
					++oceanic_collisions;
		
					plates[imap[k]].setCrust(
							x_mod,
							y_mod,
							hmap.get1D(k) - OCEANIC_BASE,
							prev_timestamp);
					
					hmap.set1D(k, hmap.get1D(k) - OCEANIC_BASE);
		
					if (hmap.get1D(k) <= 0)
					{
						imap[k] = i;
						hmap.set1D(k, this_map[j]);
						amap[k] = this_age[j];
		
						continue;
					}
				}
	
				// Record collisions to both plates. This also creates
				// continent segment at the collided location to plates.
				int this_area = plates[i].addCollision(x_mod, y_mod);
				int prev_area = plates[imap[k]].addCollision(x_mod, y_mod);
		
				// At least two plates are at same location. 
				// Move some crust from the SMALLER plate onto LARGER one.
				if (this_area < prev_area)
				{
					PlateCollision coll = new PlateCollision(
							imap[k],
							x_mod,
							y_mod,
							this_map[j] * folding_ratio);
		
					// Give some...
					hmap.set1D(k, hmap.get1D(k) + coll.crust);
					
					plates[imap[k]].setCrust(
							x_mod,
							y_mod,
							hmap.get1D(k),
							this_age[j]);
		
					// And take some.
					plates[i].setCrust(
							x_mod,
							y_mod,
							this_map[j] * (1.0f - folding_ratio),
							this_age[j]);
		
					// Add collision to the earlier plate's list.
					collisions.get(i).add(coll);
					++continental_collisions;
				}
				else
				{
					PlateCollision coll = new PlateCollision(
							i,
							x_mod,
							y_mod,
							hmap.get1D(k) * folding_ratio);
		
					plates[i].setCrust(
							x_mod,
							y_mod,
							this_map[j] + coll.crust,
							amap[k]);
		
					plates[imap[k]].setCrust(
							x_mod,
							y_mod,
							hmap.get1D(k) * (1.0f - folding_ratio),
							amap[k]);
		
					collisions.get(imap[k]).add(coll);
					++continental_collisions;
		
					// Give the location to the larger plate.
					hmap.set1D(k, this_map[j]);
					imap[k] = i;
					amap[k] = this_age[j];
				}
		    }
		}
	
	//	int total_collisions = oceanic_collisions + continental_collisions;
	//	if (total_collisions > max_collisions)
	//		max_collisions = total_collisions;
	//	printf("%5u + %5u = %5u collisions (%f %%) (max %5u (%f %%)). %c\n",
	//		oceanic_collisions, continental_collisions, total_collisions,
	//		(float)total_collisions / (float)(map_side * map_side),
	//		max_collisions, (float)max_collisions /
	//		(float)(map_side * map_side), '+' + (2 & -(iter_count & 1)));
	
		// Update the counter of iterations since last continental collision.
		last_coll_count = (last_coll_count + 1) &
			-(continental_collisions == 0 ? 1 : 0);
	
		for (int i = 0; i < num_plates; ++i)
		{
			for (PlateCollision coll : subductions.get(i))
			{
				if (DEBUG && i == coll.index)
				{
					puts("when subducting: SRC == DEST!");
					exit(1);
				}
	
				// Do not apply friction to oceanic plates.
				// This is a very cheap way to emulate slab pull.
				// Just perform subduction and on our way we go!
				plates[i].addCrustBySubduction(
					coll.wx, coll.wy, coll.crust, iter_count,
					plates[coll.index].getVelX(),
					plates[coll.index].getVelY());
			}
	
			subductions.get(i).clear();
		}
	
		for (int i = 0; i < num_plates; ++i)
		{
			for (PlateCollision coll : collisions.get(i))
			{
				int coll_count;
				float coll_ratio;
				
				CollisionInfo coll_info_i, coll_info_j;
	
				if (DEBUG && i == coll.index)
				{
					puts("when colliding: SRC == DEST!");
					exit(1);
				}
	
				// Collision causes friction. Apply it to both plates.
				plates[i].applyFriction(coll.crust);
				plates[coll.index].applyFriction(coll.crust);
	//			hmap[coll.wy * map_side + coll.wx] = 0;
	
				coll_info_i = plates[i].getCollisionInfo(coll.wx, coll.wy);
				
				coll_info_j = plates[coll.index].getCollisionInfo(coll.wx, coll.wy);
	
				// Find the minimum count of collisions between two
				// continents on different plates.
				// It's minimum because large plate will get collisions
				// from all over where as smaller plate will get just
				// a few. It's those few that matter between these two
				// plates, not what the big plate has with all the
				// other plates around it.
				coll_count = coll_info_i.count;
				coll_count -= (coll_count - coll_info_j.count) &
					-(coll_count > coll_info_j.count ? 1 : 0);
	
				// Find maximum amount of collided surface area between
				// two continents on different plates.
				// Like earlier, it's the "experience" of the smaller
				// plate that matters here.
				coll_ratio = coll_info_i.ratio;
				coll_ratio += (coll_info_j.ratio - coll_ratio) *
					(coll_info_j.ratio > coll_ratio ? 1 : 0);
	
	//			printf("min(%d, %d) = %d, max(%f, %f) = %f\n",
	//				coll_count_i, coll_count_j, coll_count,
	//				coll_ratio_i, coll_ratio_j, coll_ratio);
	
				if (coll_count > aggr_overlap_abs || coll_ratio > aggr_overlap_rel)
				{
					float amount = plates[i].aggregateCrust(
							plates[coll.index],
							coll.wx,
							coll.wy);
	
					// Calculate new direction and speed for the
					// merged plate system, that is, for the
					// receiving plate!
					plates[coll.index].collide(
							plates[i],
							coll.wx,
							coll.wy,
							amount);
				}
			}
	
			collisions.get(i).clear();
		  }
	
		// Fill divergent boundaries with new crustal material, molten magma.
		for (int y = 0, i = 0; y < (BOOL_REGENERATE_CRUST ? 1 : 0) * map_side; ++y)
		  for (int x = 0; x < map_side; ++x, ++i)
			if (imap[i] >= num_plates)
			{
				// The owner of this new crust is that neighbour plate
				// who was located at this point before plates moved.
				imap[i] = prev_imap[i];
	
				if (DEBUG && imap[i] >= num_plates)
				{
					puts("Previous index map has no owner!");
					exit(1);
				}
	
				// If this is oceanic crust then add buoyancy to it.
				// Magma that has just crystallized into oceanic crust
				// is more buoyant than that which has had a lot of
				// time to cool down and become more dense.
				amap[i] = iter_count;
				hmap.set1D(i, OCEANIC_BASE * BUOYANCY_BONUS_X);
	
				plates[imap[i]].setCrust(x, y, OCEANIC_BASE,
					iter_count);
			}
				// DEBUG!
				/*MapIndex l = new MapIndex(x, y);
				plates[imap[i]].getMapIndex(l);
				int px = (int) plates[imap[i]].left + l.x;
				int py = (int) plates[imap[i]].top + l.y;
	
				if ((py % map_side) * map_side + (px % map_side) != i)
				{
					puts("Added sea floor to odd place!");
					exit(1);
				}
			}
			else if (hmap.get1D(i) <= 0)
			{
				puts("Occupied point has no land mass!");
				exit(1);
			}*/
	
		// Add some "virginity buoyancy" to all pixels for a visual boost! :)
		for (int i = 0; i < (BUOYANCY_BONUS_X > 0 ? 1 : 0) * map_area; ++i)
		{
			// Calculate the inverted age of this piece of crust.
			// Force result to be minimum between inv. age and
			// max buoyancy bonus age.
			int crust_age = iter_count - amap[i];
			crust_age = MAX_BUOYANCY_AGE - crust_age;
			crust_age &= -(crust_age <= MAX_BUOYANCY_AGE ? 1 : 0);
	
			hmap.set1D(i, hmap.get1D(i) + (hmap.get1D(i) < CONTINENTAL_BASE ? 1 : 0) * BUOYANCY_BONUS_X *
			           OCEANIC_BASE * crust_age * MULINV_MAX_BUOYANCY_AGE);
		}
	
	/*	int i = 0;
		final int x0 = (int)plates[i].getLeft();
		final int y0 = (int)plates[i].getTop();
		final int x1 = x0 + plates[i].getWidth();
		final int y1 = y0 + plates[i].getHeight();
	
		final float*  this_map;
		final int* this_age;
		plates[i].getMap(&this_map, &this_age);
	
		// Show only plate[0]'s segments, draw everything else dark blue.
		if (iter_count < 300)
		for (int y = y0, j = 0; y < y1; ++y)
		  for (int x = x0; x < x1; ++x, ++j)
		  {
			final int x_mod = x % map_side;
			final int y_mod = y % map_side;
	
			final int k = y_mod * map_side + x_mod;
	
			if (this_map[j] < 2 * FLT_EPSILON) // No crust here...
			{
				hmap[k] = 4*FLT_EPSILON;
				continue;
			}
	
			float Q = (plates[i].segment[j] < plates[i].seg_data.size());
			hmap[k] = (this_map[j] * Q);
		  }*/
	
		//delete[] amap;
		//delete[] prev_imap;
		++iter_count;
	}
	
	@SuppressWarnings("unused")
	public void restart()
	{
		final int map_area = map_side * map_side;
		int[] amap = new int[map_area];
	
		cycle_count += max_cycles > 0 ? 1 : 0; // No increment if running for ever.
		if (cycle_count > max_cycles)
			return;
	
		// Update height map to include all recent changes.
		for (int i = 0; i < num_plates; ++i)
		{
		  final int x0 = (int)plates[i].getLeft();
		  final int y0 = (int)plates[i].getTop();
		  final int x1 = x0 + plates[i].getWidth();
		  final int y1 = y0 + plates[i].getHeight();
	
		  final float[]  this_map = plates[i].getMap();
		  final int[] this_age = plates[i].getAges();
	
		  // Copy first part of plate onto world map.
		  for (int y = y0, j = 0; y < y1; ++y)
		    for (int x = x0; x < x1; ++x, ++j)
		    {
			final int x_mod = x % map_side;
			final int y_mod = y % map_side;
	
			hmap.set(x_mod, y_mod, hmap.get(x_mod, y_mod) + this_map[j]);
			amap[y_mod * map_side + x_mod]  = this_age[j];
		    }
		}
	
		// Delete plates.
		//delete[] plates;
		plates = null;
	
		// create new plates IFF there are cycles left to run!
		// However, if max cycle count is "ETERNITY", then 0 < 0 + 1 always.
		if (cycle_count < max_cycles + (max_cycles == 0 ? 1 : 0))
		{
			//delete[] amap;
			createPlates();
			return;
		}
		//else
			//num_plates = 0;
	
		// Add some "virginity buoyancy" to all pixels for a visual boost.
		for (int i = 0; i < (BUOYANCY_BONUS_X > 0 ? 1 : 0) * map_area; ++i)
		{
			int crust_age = iter_count - amap[i];
			crust_age = MAX_BUOYANCY_AGE - crust_age;
			crust_age &= -(crust_age <= MAX_BUOYANCY_AGE ? 1 : 0);
	
			hmap.set1D(i, hmap.get1D(i) + (hmap.get1D(i) < CONTINENTAL_BASE ? 1 : 0) * BUOYANCY_BONUS_X *
			           OCEANIC_BASE * crust_age * MULINV_MAX_BUOYANCY_AGE);
		}
	
//		// This is the LAST cycle! Add some random noise to the map.
//		//int A = (map_side + 1) * (map_side + 1);
//		float[] tmp = values;
//	
//		// Shrink the fractal map by 1 pixel from right and bottom.
//		// This makes it same size as lithosphere's height map.
//		{
//			int i, j;
//			
//			for (j = 0; j < map_side; ++j)
//			for (i = 0; i < map_side; ++i)
//			{
//				hmap[j * map_side + i] = tmp[j * (map_side + 1) + i];
//			}
//		}
//	
//		float t_lowest = tmp[0], t_highest = tmp[0];
//		float h_lowest = hmap[0], h_highest = hmap[0];
//		for (int i = 1; i < map_area; ++i)
//		{
//			t_lowest = t_lowest < tmp[i] ? t_lowest : tmp[i];
//			t_highest = t_highest > tmp[i] ? t_highest : tmp[i];
//	
//			h_lowest = h_lowest < hmap[i] ? h_lowest : hmap[i];
//			h_highest = h_highest > hmap[i] ? h_highest : hmap[i];
//		}
//	
//		for (int i = 0; i < map_area; ++i)
//		{
//			// Scale to range [0, 1].
//			tmp[i] = (tmp[i] - t_lowest) / (t_highest - t_lowest);
//	
//			if (hmap[i] > CONTINENTAL_BASE)
//				hmap[i] += tmp[i] * 2;
//	//			hmap[i] = CONTINENTAL_BASE +
//	//				0.5 * (tmp[i] - 0.5) * CONTINENTAL_BASE +
//	//				0.1 * tmp[i] * (h_highest - CONTINENTAL_BASE) +
//	//				0.9 * (hmap[i] - CONTINENTAL_BASE);
//			else
//				hmap[i] = 0.8f * hmap[i] + 0.2f * tmp[i] * CONTINENTAL_BASE;
//		}
//	
//		//delete[] amap;
//		//delete[] tmp;
	}
}

