package erosion;

import src.Array2D;

public interface IErosionStyle
{
	public void setup(Array2D heightMap);
	public Point makeStartPoint();
	public float makeInitialVelocity();
	public float makeInitialSediment();
	
	public static final class Point
	{
		public int x, y;
		
		protected Point(int x, int y)
		{
			this.x = x;
			this.y = y;
		}
	}
}
