package src;

public class SubarrayDerivator implements IDerivator<Array2D>
{
	private final int x, y, width, height;
	
	public SubarrayDerivator(int x, int y, int width, int height)
	{
		this.x = x;
		this.y = y;
		this.width = width;
		this.height = height;
	}
	
	@Override
	public Array2D derive(Array2D values)
	{
		Array2D result = new Array2D(width, height);
		
		for (int i = x, newI = 0; newI < width;  ++i, ++newI)
		for (int j = y, newJ = 0; newJ < height; ++j, ++newJ)
		{
			result.set(newI, newJ, values.get(i, j));
		}
		
		return result;
	}

}
