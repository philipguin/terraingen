package squarediamond;

import src.Array2D;
import src.IPopulator;


public final class SquareDiamondArray2DPopulator implements IPopulator<Array2D>
{
	public static interface ISquareDiamondStyle
	{
		public void setup(Array2D values);
		public void onNextIteration();
		
		public float getSeedValue(int i, int j);
		public float squareValue(float topLeft, float topRight, float bottomLeft, float bottomRight);
		public float diamondValue(float left, float right, float top, float bottom);
		public float diamondValue_top(float left, float right, float bottom);
		public float diamondValue_bottom(float left, float right, float top);
		public float diamondValue_left(float right, float top, float bottom);
		public float diamondValue_right(float left, float top, float bottom);
	}
	
	private final ISquareDiamondStyle style;
	
	public SquareDiamondArray2DPopulator(ISquareDiamondStyle style)
	{
		this.style = style;
	}
	
	@Override
	public final void populate(Array2D values)
	{
		if (values.width != values.height)
			throw new Error("values' width was not equal to height!");
		
		int test = values.width - 2;
		
		while (test > 0)
		{
			if ((test & 1) == 0)
				throw new Error("values was not a power of 2 plus 1!");
			
			test >>= 1;
		}
		
		int i, j;
		int width = values.width;
		int jump = width - 1;
		int halfJump = jump >> 1;

		for (j = 0; j < width; j += jump)
		for (i = 0; i < width; i += jump)
		{
			values.set(i, j, style.getSeedValue(i, j));
		}
		
		style.setup(values);
		
		while (jump > 1)
		{
			style.onNextIteration();
			
			//Ref Square
			for (j = halfJump; j < width; j += jump)
			for (i = halfJump; i < width; i += jump)
			{
				values.set(i, j, style.squareValue(
						values.get(i - halfJump, j - halfJump),
						values.get(i + halfJump, j - halfJump),
						values.get(i - halfJump, j + halfJump),
						values.get(i + halfJump, j + halfJump)));
			}
			
			//Ref Diamond
			for (j = 0; j < width; j += halfJump)
			{
				for (i = halfJump; i < width; i += jump)
				{
					if (j == 0)
					{
						values.set(i, j, style.diamondValue_top(
								values.get(i - halfJump, j),
								values.get(i + halfJump, j),
								values.get(i, j + halfJump)));
					}
					else if (j == width-1)
					{
						values.set(i, j, style.diamondValue_bottom(
								values.get(i - halfJump, j),
								values.get(i + halfJump, j),
								values.get(i, j - halfJump)));
					}
					else
					{
						values.set(i, j, style.diamondValue(
								values.get(i - halfJump, j),
								values.get(i + halfJump, j),
								values.get(i, j - halfJump),
								values.get(i, j + halfJump)));
					}
				}
				
				j += halfJump;
				
				if (j >= width)
					break;
				
				for (i = 0; i < width; i += jump)
				{
					if (i == 0)
					{
						values.set(i, j, style.diamondValue_left(
								values.get(i + halfJump, j),
								values.get(i, j - halfJump),
								values.get(i, j + halfJump)));
					}
					else if (i == width - 1)
					{
						values.set(i, j, style.diamondValue_right(
								values.get(i - halfJump, j),
								values.get(i, j - halfJump),
								values.get(i, j + halfJump)));
					}
					else
					{
						values.set(i, j, style.diamondValue(
								values.get(i - halfJump, j),
								values.get(i + halfJump, j),
								values.get(i, j - halfJump),
								values.get(i, j + halfJump)));
					}
				}
			}

			jump = halfJump;
			halfJump >>= 1;
		}
	}

}
