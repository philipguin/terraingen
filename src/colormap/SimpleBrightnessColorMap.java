package colormap;
import java.awt.Color;



public class SimpleBrightnessColorMap implements Mapper<Float, Color>
{
	private final float hue, saturation;
	
	public SimpleBrightnessColorMap(float hue, float saturation)
	{
		this.hue = hue;
		this.saturation = saturation;
	}
	
	@Override
	public Color map(Float value)
	{
		return new Color(Color.HSBtoRGB(hue, saturation, value));
	}
}
